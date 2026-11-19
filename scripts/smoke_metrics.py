"""Fresh all-locus technical metric integration; no scientific calibration.

Separate matched tryptic/chymotryptic acquisitions are pooled after search.
Observed native ARDA sequences and official constant references alone form the DB.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.smoke_reads import selected_modes, expected_pair_count


def preflight(root, cpus):
    if root.exists():
        raise FileExistsError('output already exists; choose a fresh directory')
    if cpus < 1:
        raise ValueError('CPU budget must be positive')


def aggregate_coverage(frame, keys):
    """Keep absent/undefined region populations distinct from measured zero."""
    import polars as pl
    result = frame.group_by(keys).agg(
        pl.len().alias('represented_chains'),
        pl.col('coverage').count().alias('defined_chains'),
        pl.col('region_status').first().alias('region_status'),
        (pl.col('covered_residues')>0).sum().alias('covered_chains'),
        pl.col('region_length').sum().alias('residue_denominator'),
        pl.col('covered_residues').sum().alias('covered_residues'),
        pl.col('coverage').mean().alias('macro_chain_coverage'))
    return result.with_columns(
        *[pl.when(pl.col('defined_chains')>0).then(pl.col(c)).otherwise(None).alias(c)
          for c in ('covered_chains','residue_denominator','covered_residues')],
        pl.when(pl.col('defined_chains')>0).then(pl.lit('defined')).otherwise(pl.col('region_status')).alias('status')
    ).with_columns(pl.when(pl.col('residue_denominator')>0)
        .then(pl.col('covered_residues')/pl.col('residue_denominator'))
        .otherwise(None).alias('pooled_residue_coverage'))


def aggregate_gene_coverage(coverage, metadata):
    """Preserve native bins; attach scalar V labels without duplicating aliases."""
    import polars as pl
    fields=['canonical_v_gene','canonical_v_match_status','canonical_v_gene_count']
    lookup=metadata.select('locus',pl.col('v_call').alias('gene_call'),*fields).unique()
    genes=aggregate_coverage(coverage.filter(pl.col('region').is_in(['V','D','J','C'])),
                             ['locus','region','gene_call','evidence_level'])
    return genes.join(lookup,on=['locus','gene_call'],how='left',validate='m:1').with_columns(
        pl.when(pl.col('region')=='V').then(pl.col('canonical_v_match_status').fill_null('not_evaluated'))
          .otherwise(pl.lit('not_applicable')).alias('canonical_v_match_status'),
        *[pl.when(pl.col('region')=='V').then(pl.col(f)).otherwise(None).alias(f)
          for f in ['canonical_v_gene','canonical_v_gene_count']]
    ).sort(['locus','region','gene_call','evidence_level'])


def component_counts(metadata):
    """Count observed physical rows; latent ancestors never enter these denominators."""
    import polars as pl
    return metadata.group_by(['component','locus']).agg(
        pl.len().alias('observed_chain_count'),
        pl.col('antibody_id').n_unique().alias('observed_antibody_count'),
        pl.col('family_id').n_unique().alias('observed_family_count')
    ).sort(['component','locus'])


def component_recovery(identity):
    """Stratify global identity outcomes without resolving cross-component aliases."""
    import polars as pl
    return identity.group_by('component').agg(
        pl.len().alias('observed_heavy_rows'),
        pl.col('variable_mass_class_id').n_unique().alias('variable_mass_classes'),
        pl.col('variable_mass_class_id').filter(pl.col('full_variable_compatible_candidate_count')>0)
          .n_unique().alias('observable_classes'),
        pl.col('variable_mass_class_id').filter(pl.col('class_recovered_at20'))
          .n_unique().alias('recovered_classes_at20')
    ).with_columns(
        (pl.col('recovered_classes_at20')/pl.col('variable_mass_classes')).alias('class_recovery_at20'),
        pl.when(pl.col('observable_classes')>0)
          .then(pl.col('recovered_classes_at20')/pl.col('observable_classes'))
          .otherwise(None).alias('conditional_class_recovery_at20')
    ).sort('component')


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--geometry', choices=('both','amplicon','rnaseq'))
    p.add_argument('--cpus', type=int)
    p.add_argument('--reference-dir', type=Path, help='checksum-pinned official raw records; never generated outputs')
    population = p.add_mutually_exclusive_group()
    population.add_argument('--families', type=int, help='positive even family count; activates population generator')
    population.add_argument('--total-antibodies', type=int, help='observed paired records, excluding latent ancestors')
    p.add_argument('--memory-fraction', type=float, choices=(.25,.5,.75))
    p.add_argument('--heavy-constant-gene', choices=('IGHG1','IGHM','IGHG3','IGHA1','IGHE'))
    p.add_argument('--additional-reference-dir', type=Path, help='pinned raw records for nondefault heavy constants')
    p.add_argument('--transcript-reference-dir', type=Path, help='enable separately versioned leader-imputed RNA model')
    p.add_argument('--pairs-per-chain', type=int)
    p.add_argument('--mutation-rate', type=float)
    p.add_argument('--protein-concentration', choices=('equal','head80'))
    p.add_argument('--seed', type=int)
    p.add_argument('--read-seed', type=int, default=2468)
    p.add_argument('--gradient-seconds', type=float, default=120.)
    p.add_argument('--ms2-per-bin', type=int, default=5)
    p.add_argument('--purpose', choices=('technical','resource-profile','mechanistic'), default='technical')
    p.add_argument('--mssimulator',type=Path,help='explicit MSSimulator executable override')
    p.add_argument('--stage-timeout',type=float)
    args = p.parse_args(argv)
    if not 0 <= args.read_seed < 2**32:p.error('--read-seed must be in [0,2**32)')
    if args.transcript_reference_dir is not None:
        if (args.geometry != 'rnaseq' or args.total_antibodies is None or
                not args.transcript_reference_dir.is_dir()):
            p.error('transcript model requires --geometry rnaseq, --total-antibodies and existing raw transcript-reference-dir')
    args.stage_timeout_explicit = args.stage_timeout is not None
    if args.stage_timeout is None: args.stage_timeout = 300.
    if not math.isfinite(args.stage_timeout) or args.stage_timeout<=0:
        p.error('--stage-timeout must be finite and positive')
    if args.mssimulator is not None:
        args.mssimulator=args.mssimulator.resolve()
        if not args.mssimulator.is_file() or not os.access(args.mssimulator,os.X_OK):
            p.error('--mssimulator must be an executable file')
    if args.families is None and args.total_antibodies is None:
        if any(getattr(args, key) is not None for key in ('pairs_per_chain','mutation_rate','protein_concentration','seed','heavy_constant_gene','additional_reference_dir','memory_fraction')):
            p.error('population settings require --families or --total-antibodies')
    else:
        if args.families is not None and (not 0 < args.families < 10**10 or args.families % 2):
            p.error('--families must be positive, even and below 10000000000')
        if args.total_antibodies is not None:
            if args.total_antibodies <= 0 or args.total_antibodies >= 3*10**10 or args.memory_fraction is None:
                p.error('mixture requires positive --total-antibodies and --memory-fraction')
            memory = args.total_antibodies * args.memory_fraction
            if memory % 6 or (args.total_antibodies-memory) % 2:
                p.error('mixture requires whole balanced memory families and naive singleton pairs')
        elif args.memory_fraction is not None:
            p.error('--memory-fraction requires --total-antibodies')
        for key, default in [('pairs_per_chain',100),('mutation_rate',.01),('protein_concentration','equal'),('seed',20261110),('heavy_constant_gene','IGHG1')]:
            if getattr(args, key) is None:
                setattr(args, key, default)
        if args.pairs_per_chain <= 0:
            p.error('--pairs-per-chain must be positive')
        if not math.isfinite(args.mutation_rate) or not 0 <= args.mutation_rate <= 1:
            p.error('--mutation-rate must be finite in [0,1]')
        if not 0 <= args.seed < 10**14:
            p.error('--seed must be in [0,100000000000000)')
        if args.total_antibodies is not None:
            if args.heavy_constant_gene != 'IGHG1' or args.additional_reference_dir is None:
                p.error('mixture requires default memory IGHG1 and --additional-reference-dir for naive IGHM')
        elif (args.heavy_constant_gene != 'IGHG1') != (args.additional_reference_dir is not None):
            p.error('nondefault heavy constants require --additional-reference-dir; IGHG1 uses the base references')
        if args.additional_reference_dir is not None and not args.additional_reference_dir.is_dir():
            p.error('--additional-reference-dir must be an existing raw reference directory')
    if not math.isfinite(args.gradient_seconds) or args.gradient_seconds <= 0:
        p.error('--gradient-seconds must be finite and positive')
    if args.ms2_per_bin < 1:
        p.error('--ms2-per-bin must be positive')
    return args


def generate_inputs(args, repo, root, run):
    from smoke_all import prepare_references
    root.parent.mkdir(parents=True, exist_ok=True)
    generation = repo/'.tool-envs/generation/bin/python'
    if args.families is None and args.total_antibodies is None:
        prepare_references(repo, root, all_loci=True, source_dir=args.reference_dir)
        run('generation', generation, 'smoke_generate.py', '--output', root, '--all-loci')
    else:
        def generate_population(reference_dir):
            population = (['--total-antibodies', args.total_antibodies, '--memory-fraction', args.memory_fraction]
                          if args.total_antibodies is not None else ['--families', args.families])
            additional = ['--additional-reference-dir', args.additional_reference_dir.resolve()] if args.additional_reference_dir else []
            transcript = (['--transcript-reference-dir', args.transcript_reference_dir.resolve(),
                           '--transcript-threads', __import__('vdjaqp').resolve_cpu_budget(args.cpus)]
                          if args.transcript_reference_dir else [])
            run('generation', generation, 'generate_population.py', '--output', root,
                '--reference-dir', reference_dir, *population,
                '--pairs-per-chain', args.pairs_per_chain, '--mutation-rate', args.mutation_rate,
                '--protein-concentration', args.protein_concentration, '--seed', args.seed,
                '--heavy-constant-gene', args.heavy_constant_gene, *additional, *transcript)
        if args.reference_dir:
            generate_population(args.reference_dir.resolve())
        else:
            with TemporaryDirectory(prefix='vdjaqp-raw-references-') as temporary:
                raw_root = Path(temporary)
                prepare_references(repo, raw_root, all_loci=True)
                generate_population(raw_root/'references')


def constant_search_fasta(root, *, population):
    """Use declared reference proteins, with no access to generated variable truth."""
    if population:
        return (root/'constant_reference_proteins.fasta').read_text()
    return ''.join('>constant_'+accession+'\n'+''.join(
        (root/'references'/f'{accession}.fasta').read_text().splitlines()[1:])+'\n'
        for accession in ('P01857-1','P01834','P0DOY2'))

def load_observed(path, geometry):
    """Read only a completed native ARDA export, checking its final receipt."""
    import polars as pl
    receipt=json.loads((path/'receipt.json').read_text())
    if receipt.get('status')!='complete':
        raise ValueError('native ARDA receipt must be complete')
    if receipt['source_format']!='arda_fastq_native' or receipt['geometry']!=geometry:
        raise ValueError('native ARDA receipt source/geometry mismatch')
    source_pins=receipt.get('source_input_sha256')
    if not isinstance(source_pins,dict) or not source_pins:
        raise ValueError('native ARDA source pins are required')
    verified=[];unavailable=[]
    for name,pin in source_pins.items():
        try:
            with Path(name).open('rb') as handle:
                digest=hashlib.file_digest(handle,'sha256').hexdigest()
        except (FileNotFoundError,PermissionError):
            unavailable.append(name)
            continue
        if digest!=pin:raise ValueError('native ARDA source hash mismatch: '+name)
        verified.append(name)
    frames={'source_verification':dict(verified_local_sources=verified,
        unavailable_sources_not_verified=unavailable)}
    for name in ('catalog','candidates','search_entries'):
        artifact=path/(name+'.parquet')
        with artifact.open('rb') as handle:
            digest=hashlib.file_digest(handle,'sha256').hexdigest()
        if digest!=receipt['artifact_sha256'][artifact.name]:
            raise ValueError('native ARDA artifact hash mismatch: '+name)
        frame=pl.read_parquet(artifact)
        if frame.height!=receipt['counts'][name]:
            raise ValueError('native ARDA artifact count mismatch: '+name)
        if (frame.schema.get('background_member')!=pl.Boolean or
            frame['background_member'].null_count() or not frame['background_member'].all()):
            raise ValueError('native ARDA background membership must be explicit and true')
        frames[name]=frame
    return frames


def main():
    args = parse_args()
    modes = selected_modes(args.geometry or 'both')
    geometry_flags = ('--geometry', args.geometry) if args.geometry is not None else ()
    repo = Path(__file__).resolve().parents[1]
    import vdjaqp
    cpus = vdjaqp.resolve_cpu_budget(args.cpus)
    root = args.output.resolve()
    preflight(root, cpus)
    import numpy as np
    import polars as pl
    from arda.annotate.contig import reannotate_contigs
    from vdjaqp.simulation_truth import regions_from_annotations
    from vdjaqp.simulation_metrics import exact_peptide_placements, region_coverage, psm_correctness, identity_metrics
    from vdjaqp.readout import read_idxml
    from vdjaqp.pipeline import Evidence
    from vdjaqp.errormodel import ErrorModel, mass_prior, _E_DB_DEFAULT
    from vdjaqp.bench import rank_targets
    started = time.perf_counter()
    stages = []
    def run(label, python, script, *arguments):
        command = [str(python), str(repo/'scripts'/script), *map(str, arguments)]
        if script in ('smoke_reads.py','smoke_repertoire.py','smoke_proteomics.py') and args.stage_timeout_explicit:
            command.extend(['--stage-timeout',str(args.stage_timeout)])
        if script=='smoke_proteomics.py':
            if args.mssimulator is not None:command.extend(['--mssimulator',str(args.mssimulator)])
        begin = time.perf_counter()
        print('Starting '+label, flush=True)
        # Generation requires an output containing raw references only.
        logpath = root.parent/(root.name+'-'+label+'.log') if label == 'generation' else root/(label+'.log')
        with logpath.open('w') as log:
            subprocess.run(command, cwd=repo, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=2*args.stage_timeout+60 if script=='smoke_proteomics.py' else args.stage_timeout)
        if label == 'generation': logpath.rename(root/'generation.log')
        stages.append(dict(stage=label, command=command, seconds=time.perf_counter()-begin))
        print('Completed '+label, flush=True)
    def save(frame, name):
        frame.write_parquet(root/(name+'.parquet'), compression='zstd')
    generate_inputs(args, repo, root, run)
    generation = repo/'.tool-envs/generation/bin/python'
    expected_pairs = expected_pair_count(root/'read_counts.tsv')
    truth = json.loads((root/'truth.json').read_text())
    assert expected_pairs == sum(c['rna_pairs'] for c in truth['chains'])
    if args.total_antibodies is not None:
        assert len(truth['chains']) == 2*args.total_antibodies
        assert len({c['antibody_id'] for c in truth['chains']}) == args.total_antibodies
        assert expected_pairs == 2*args.total_antibodies*args.pairs_per_chain
        heavy = [c for c in truth['chains'] if c['locus']=='IGH']
        assert sum(c['component']=='memory' for c in heavy) == args.total_antibodies*args.memory_fraction
        assert sum(c['component']=='naive' for c in heavy) == args.total_antibodies*(1-args.memory_fraction)
    elif args.families is not None:
        assert len(truth['chains']) == 6*args.families
        assert expected_pairs == 6*args.families*args.pairs_per_chain
    run('reads', generation, 'smoke_reads.py', '--root', root, '--cpus', cpus, '--read-seed', args.read_seed, *geometry_flags)
    read_runs = json.loads((root/'read_run_manifest.json').read_text())
    assert {r['scenario'] for r in read_runs} == set(modes) and len(read_runs) == len(modes)
    assert all(r['pairs'] == expected_pairs for r in read_runs)
    for mode in modes:
        run(mode+'_repertoire', sys.executable, 'smoke_repertoire.py', mode, '--root', root, '--threads', cpus)
    chains = json.loads((root/'truth.json').read_text())['chains']
    annotation_started = time.perf_counter()
    annotations = reannotate_contigs([(c['chain_id'], c['variable_nt']) for c in chains], seqtype='nt', threads=cpus)
    stages.append(dict(stage='evaluation_truth_annotation',seconds=time.perf_counter()-annotation_started,input_chains=len(chains),threads=cpus,seqtype='nt'))
    (root/'truth_annotations.json').write_text(json.dumps(annotations)+'\n')
    metadata, regions = regions_from_annotations(chains, annotations)
    save(metadata, 'chain_truth'); save(regions, 'truth_regions')
    metadata.write_csv(root/'chain_truth.tsv', separator='\t')
    repertoire = metadata.with_columns((pl.col('rna_pairs')/pl.col('rna_pairs').sum()).alias('rna_pool_fraction'),(pl.col('rna_pairs')*pl.col('secretion_factor')).alias('serum_weight')).with_columns((pl.col('serum_weight')/pl.col('serum_weight').sum()).alias('serum_pool_fraction'),pl.len().over('family_id').alias('family_chain_count'),pl.col('antibody_id').n_unique().over('family_id').alias('family_antibody_count'),pl.lit(len(chains)).alias('physical_chain_denominator'))
    save(repertoire, 'repertoire_metrics')
    if args.total_antibodies is not None:
        save(component_counts(metadata), 'component_counts')
    observed = {mode: load_observed(root/(mode+'_observed'), mode) for mode in modes}
    catalog, candidates, entries = (pl.concat([observed[mode][key] for mode in modes], how='diagonal_relaxed')
                                    for key in ('catalog','candidates','search_entries'))
    save(catalog, 'observation_catalog'); save(candidates, 'candidates'); save(entries, 'search_entries')
    physical = metadata.select('chain_id','protein_aa')
    levels = ['native_ms2_contributor_compatibility','identified_peptide_compatibility','arda_observed_exact_compatibility']
    summaries = {}
    discarded = {}
    constant_fasta = constant_search_fasta(root, population=args.families is not None or args.total_antibodies is not None)
    for mode in modes:
        # Geometry membership comes from native ARDA exports, never simulation truth.
        obs = observed[mode]['search_entries']
        clones = observed[mode]['candidates']
        database = root/(mode+'_search.fasta')
        text = ''.join('>'+r['entry_id']+'\n'+r['seq']+'\n' for r in obs.to_dicts())
        text += constant_fasta
        database.write_text(text)
        save(clones, mode+'_candidates')
        calls, truth_links, native_queries = [], [], []
        spectra = []
        for enzyme in ('Trypsin','Chymotrypsin'):
            label = mode+'_'+enzyme
            destination = root/label
            run(label, sys.executable, 'smoke_proteomics.py', root/'proteins.fasta', database, destination, '--enzyme', enzyme, '--threads', cpus, '--keep-spectra', '--gradient-seconds', args.gradient_seconds, '--ms2-per-bin', args.ms2_per_bin)
            native = json.loads((destination/'ms2_truth.json').read_text())
            peptides = read_idxml(destination/'psms.idXML', run_id=label).with_columns(pl.lit(enzyme).alias('enzyme'))
            per, summary = psm_correctness(peptides, native)
            assert per.filter(pl.col('status')!='evaluated').height == 0
            assert set(peptides['native_spectrum_reference']) <= {r['native_spectrum_reference'] for r in native}
            save(per, label+'_psm_correctness'); save(summary, label+'_readout_metrics')
            calls.append(peptides)
            spectra.append(destination/'spectra.mzML')
            for event in native:
                for contributor in event['contributors']:
                    seq = contributor['peptide_aa']
                    native_queries.append(seq)
                    for chain in contributor['protein_accessions']:
                        truth_links.append({'seq':seq, 'chain_id':chain})
        peptides = pl.concat(calls, how='vertical_relaxed')
        save(peptides, mode+'_peptides')
        universe = sorted(set(native_queries+peptides['seq'].to_list()+obs['seq'].to_list()))
        queries = pl.DataFrame({'peptide_id':[f'q{i}' for i in range(len(universe))], 'seq':universe},schema={'peptide_id':pl.String,'seq':pl.String})
        short = queries.filter(pl.col('seq').str.len_chars()<6)
        save(short, mode+'_unsearchable_queries')
        queries = queries.filter(pl.col('seq').str.len_chars()>=6)
        hits = exact_peptide_placements(physical, queries, threads=cpus).join(queries,on='peptide_id')
        save(hits, mode+'_exact_placements')
        links = pl.DataFrame(truth_links, schema={'seq':pl.String,'chain_id':pl.String}).unique()
        evidence = [hits.join(links,on=['seq','chain_id'],how='semi').with_columns(pl.lit(levels[0]).alias('evidence_level'))]
        for sequences, level in ((peptides['seq'],levels[1]),(obs['seq'],levels[2])):
            evidence.append(hits.filter(pl.col('seq').is_in(sequences.implode())).with_columns(pl.lit(level).alias('evidence_level')))
        coverage = region_coverage(regions,pl.concat(evidence).select('chain_id','start','end','evidence_level'),levels)
        save(coverage, mode+'_chain_region_coverage')
        coverage.write_csv(root/(mode+'_chain_region_coverage.tsv'), separator='\t')
        genes = aggregate_gene_coverage(coverage, metadata)
        save(genes,mode+'_gene_coverage')
        genes.write_csv(root/(mode+'_gene_coverage.tsv'), separator='\t')
        save(aggregate_coverage(coverage,['locus','region','evidence_level']),mode+'_locus_coverage')
        if args.total_antibodies is not None:
            component_coverage = coverage.join(metadata.select('chain_id','component'),on='chain_id',how='left',validate='m:1')
            save(aggregate_coverage(component_coverage,['component','locus','region','evidence_level']),mode+'_component_coverage')
        if clones.height:
            model = ErrorModel(np.array([0.,1.]),np.array([.1,.01]),mass_prior(),_E_DB_DEFAULT)
            ev = Evidence.from_frames(clones,peptides,error_model=model)
            scores = ev.score()
            assert scores.height == clones.height
            assert np.isfinite(scores['score'].to_numpy()).all()
            alignments = ev.matches.height
            del ev
        else:
            scores = pl.DataFrame(schema={'clone_id':pl.String,'score':pl.Float64,'confirmed_nats':pl.Float64})
            alignments = 0
        save(scores,mode+'_scores')
        ranks = rank_targets(scores,'score',set(scores['clone_id']))
        save(ranks,mode+'_ranks')
        association = clones.select('clone_id','seq').join(hits.filter(pl.col('chain_id').is_in(metadata.filter(pl.col('locus')=='IGH')['chain_id'].implode())),on='seq').select('clone_id','chain_id','start','end').unique().join(ranks,on='clone_id')
        save(association,mode+'_candidate_associations')
        identity, candidate_identity, recovery = identity_metrics(metadata, clones, ranks, association)
        save(identity,mode+'_identity_population')
        save(candidate_identity,mode+'_candidate_identity')
        save(recovery,mode+'_identity_metrics')
        if args.total_antibodies is not None:
            save(component_recovery(identity),mode+'_component_recovery')
        summaries[mode] = dict(candidates=clones.height,search_entries=obs.height,psm_assignments=peptides.height,exact_placements=hits.height,coverage_rows=coverage.height,association_rows=association.height,alignment_rows=alignments,true_heavy_denominator=metadata.filter(pl.col('locus')=='IGH').height,identity_status='full-variable I/L-equivalent mass classes; existing >=20 confirmed-nat criterion',identity_metrics=recovery.to_dicts())
        if args.total_antibodies is not None:
            summaries[mode]['component_identity_metrics'] = component_recovery(identity).to_dicts()
        for path in spectra:
            discarded[str(path.relative_to(root))] = dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bytes=path.stat().st_size)
            path.unlink()
    files = [f for f in root.rglob('*') if f.is_file()]
    size = sum(f.stat().st_size for f in files)
    assert size < 10_000_000_000
    scenario_inputs = dict(purpose=args.purpose, families=args.families if args.total_antibodies is not None or args.families is not None else 2,
        geometry=args.geometry or 'both',
        total_antibodies=args.total_antibodies, memory_fraction=args.memory_fraction,
        antibodies=len({c['antibody_id'] for c in chains}), chains=len(chains),
        expected_read_pairs_per_geometry=expected_pairs, actual_read_pairs={r['scenario']:r['pairs'] for r in read_runs},
        seed=truth['seed'], read_seed=args.read_seed, mutation_rate=args.mutation_rate if args.families is not None or args.total_antibodies is not None else .01,
        protein_concentration=args.protein_concentration if args.families is not None or args.total_antibodies is not None else 'technical fixture weights',
        heavy_constant_gene=args.heavy_constant_gene or 'IGHG1',
        naive_constant_gene='IGHM' if args.total_antibodies is not None else None,
        pairs_per_chain=args.pairs_per_chain, gradient_seconds=args.gradient_seconds, ms2_per_bin=args.ms2_per_bin,
        reconstruction_claim=False,stage_timeout_seconds=args.stage_timeout,
        mssimulator_override=str(args.mssimulator) if args.mssimulator else None)
    source_files = [*repo.glob('vdjaqp/*.py'), *repo.glob('scripts/smoke*.py'), repo/'scripts/generate_population.py',
                    repo/'scripts/genairr_gene_map.py', repo/'scripts/validate_constant_candidates.py',
                    *[f for f in (repo/'envs').rglob('*') if f.is_file()]]
    if args.transcript_reference_dir is not None:
        scenario_inputs['rna_model'] = 'leader_imputed_v1'
        source_files.extend([repo/'scripts/build_transcript_context.py', repo/'scripts/audit_transcript_context.py'])
    report = dict(observed_source_verification={mode:observed[mode]['source_verification'] for mode in modes},scenario_inputs=scenario_inputs,raw_reference_source=str(args.reference_dir.resolve()) if args.reference_dir else 'official manifests download',cpus=cpus,seconds=time.perf_counter()-started,stages=stages,scenarios=summaries,retained_bytes=size,discarded_spectra=discarded,artifact_sha256={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},software_sha256={str(f.relative_to(repo)):hashlib.sha256(f.read_bytes()).hexdigest() for f in source_files},peak_child_rss_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
    (root/'metric_integration.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(summaries,indent=2))


if __name__ == '__main__':
    main()
