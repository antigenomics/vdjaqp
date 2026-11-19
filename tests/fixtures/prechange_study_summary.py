"""Consolidate the twelve approved completed builds' final tables for R figures.

No generation, annotation, sequence search, scoring, fitting or significance tests.
Input bundles and output tables remain private requested deliverables.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from vdjaqp import resolve_cpu_budget

CASES = {
    'M0': (.5, .01, 'equal', 5), 'M1': (.5, .01, 'equal', 10),
    'M2': (.5, .01, 'head80', 5), 'M3': (.5, .03, 'equal', 5),
    'M4': (.25, .01, 'equal', 5), 'M5': (.75, .01, 'equal', 5),
}
GEOMETRIES = ('amplicon', 'rnaseq')
MONITOR_PINS = {
    '57d57da451927c7f695fdb2f1d0387ac6a7fe014427bdb9281b36c07919b040a',
    '39895acf682ef491f0d01ff178a83fa73ed3dbdef458ec229f1e63300d147acf',
}
SCIENTIFIC_ENTRYPOINTS = {'scripts/smoke_metrics.py', 'scripts/smoke_reads.py',
    'scripts/smoke_assemble.py', 'scripts/smoke_proteomics.py', 'vdjaqp/simulation_truth.py'}
EXISTING = ('chain_region_coverage', 'locus_coverage', 'component_coverage',
            'gene_coverage', 'identity_population', 'identity_metrics', 'component_recovery')


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def inside(root, name):
    p = Path(name)
    if p.is_absolute() or '..' in p.parts or not p.parts:
        raise ValueError('bundle paths must be relative and confined')
    target = root/p
    if not target.resolve().is_relative_to(root) or any(x.is_symlink() for x in [target, *target.parents] if x.is_relative_to(root)):
        raise ValueError('bundle path escapes or traverses symlink')
    return target


def inventory(root):
    expected = {f'{case}-{mode}' for case in CASES for mode in GEOMETRIES}
    if (not root.is_dir() or {p.name for p in root.iterdir()} != expected
            or any(not p.is_dir() or p.is_symlink() for p in root.iterdir())):
        raise ValueError('input root must contain exactly the twelve approved case directories')
    return [(case, mode, root/f'{case}-{mode}') for case in CASES for mode in GEOMETRIES]


def validate_bundle(root, case, mode):
    manifest = json.loads(inside(root, 'retained_table_manifest.json').read_text())
    records = manifest['artifacts']
    mapped = {r['file']:r for r in records}
    if len(mapped) != len(records): raise ValueError('duplicate manifest artifact')
    required = {'metric_integration.json', 'resource.json', 'entrypoints.sha256',
                'repertoire_metrics.parquet', 'truth_regions.parquet', 'component_counts.parquet',
                'contig_catalog.parquet', 'assembly/arda_contigs.json.gz',
                mode+'_candidate_identity.parquet'}
    required |= {mode+'_'+name+'.parquet' for name in EXISTING}
    required |= {mode+'_'+enzyme+'_'+name+'.parquet' for enzyme in ('Trypsin','Chymotrypsin')
                 for name in ('psm_correctness','readout_metrics')}
    if not required <= mapped.keys(): raise ValueError('required final tables missing from manifest')
    report = json.loads(inside(root, 'metric_integration.json').read_text())
    for name, record in mapped.items():
        path = inside(root, name)
        if path.stat().st_size != record['bytes'] or digest(path) != record['sha256']:
            raise ValueError('retained artifact hash/size mismatch: '+name)
        if 'compression' in record:
            if not name.endswith('.gz'): raise ValueError('unsupported retained compression')
            with gzip.open(path, 'rb') as stream: raw = stream.read()
            if len(raw) != record['source_bytes'] or hashlib.sha256(raw).hexdigest() != record['source_sha256']:
                raise ValueError('compressed source hash mismatch')
            source_hash = record['source_sha256']
        else:
            source_hash = record['sha256']
        source = record['source_path']
        inside(root, source)
        if source in report['artifact_sha256'] and report['artifact_sha256'][source] != source_hash:
            raise ValueError('retained source differs from producer report')
        if name not in ('metric_integration.json','resource.json','entrypoints.sha256') and source not in report['artifact_sha256']:
            raise ValueError('final table is not listed by producer')
    params = report['scenario_inputs']
    fraction, mutation, concentration, capacity = CASES[case]
    expected = dict(geometry=mode, memory_fraction=fraction, mutation_rate=mutation,
                    protein_concentration=concentration, ms2_per_bin=capacity, seed=123,
                    pairs_per_chain=100, gradient_seconds=1200, purpose='mechanistic',
                    reconstruction_claim=False, antibodies=1200, chains=2400,
                    total_antibodies=1200, actual_read_pairs={mode:240000})
    if any(params.get(k) != v for k,v in expected.items()):
        raise ValueError('case parameters differ from approved inventory')
    if set(report['scenarios']) != {mode}: raise ValueError('mixed RNA geometries')
    resources = json.loads((root/'resource.json').read_text())
    if (resources['failure'] is not None or resources['exit_code'] != 0
            or resources['cpu_budget'] != 16 or resources['live_thread_limit'] != 64
            or resources['peak_live_threads'] > 64
            or resources['peak_sampled_summed_rss_bytes'] > 24_000_000_000
            or resources['wall_seconds'] >= 1800
            or resources['peak_sampled_project_allocated_bytes'] > 40_000_000_000):
        raise ValueError('case did not pass its resource contract')
    expected_spectra = {f'{mode}_{e}/spectra.mzML' for e in ('Trypsin','Chymotrypsin')}
    if set(report['discarded_spectra']) != expected_spectra:
        raise ValueError('producer lacks complete two-protease spectra cleanup record')
    # The monitor is recorded by the job, outside the scientific producer's hash list.
    entrypoints = {}
    for line in (root/'entrypoints.sha256').read_text().splitlines():
        sha, name = line.split(maxsplit=1)
        name = name.lstrip('*')
        if name in entrypoints: raise ValueError('duplicate entrypoint hash')
        entrypoints[name] = sha
        if name == 'scripts/profile_build.py':
            if sha not in MONITOR_PINS: raise ValueError('unknown resource-monitor source hash')
            if name in report['software_sha256'] and report['software_sha256'][name] != sha:
                raise ValueError('monitor hash conflicts with producer report')
        elif report['software_sha256'].get(name) != sha:
            raise ValueError('entrypoint hash differs from producer report')
    if not SCIENTIFIC_ENTRYPOINTS <= entrypoints.keys() or 'scripts/profile_build.py' not in entrypoints:
        raise ValueError('required scientific or monitoring entrypoint missing')
    return mapped, report, {'scripts/profile_build.py':entrypoints['scripts/profile_build.py']}


def gene_frequencies(frame, fields, population, namespace):
    import polars as pl
    tables = []
    for segment, column in fields.items():
        part = frame.select('locus', pl.col(column).alias('gene_call'))
        if segment=='V' and 'canonical_v_gene' in frame.columns:
            part = part.with_columns(frame['canonical_v_gene'], frame['canonical_v_match_status'],
                frame['canonical_v_gene_count'] if 'canonical_v_gene_count' in frame.columns else
                pl.lit(None,dtype=pl.Int64).alias('canonical_v_gene_count'))
        else:
            part = part.with_columns(pl.lit(None,dtype=pl.String).alias('canonical_v_gene'),
                                     pl.lit('not_evaluated' if segment=='V' else 'not_applicable').alias('canonical_v_match_status'),
                                     pl.lit(None,dtype=pl.Int64).alias('canonical_v_gene_count'))
        part = part.with_columns((pl.col('gene_call').is_null() | (pl.col('gene_call')=='')).alias('missing_call'))
        totals = part.group_by('locus').len(name='denominator')
        counts = part.group_by(part.columns).len(name='numerator').join(totals,on='locus',nulls_equal=True)
        tables.append(counts.with_columns(pl.lit(segment).alias('segment'),pl.lit(population).alias('population'),
            pl.lit(namespace).alias('gene_namespace'),pl.lit('record_fraction').alias('unit'),
            pl.col('denominator').alias('n_observations'),
            (pl.col('numerator')/pl.col('denominator')).alias('fraction')))
    return pl.concat(tables,how='vertical_relaxed')


def summarize_tables(root, case, mode, records):
    import polars as pl
    from vdjaqp.airr_ms_metrics import eligible_species
    output = {}
    def tag(frame, source, population):
        if any(dtype.is_nested() for dtype in frame.schema.values()): raise ValueError('nested figure table')
        return frame.with_columns(pl.lit(case).alias('case'),pl.lit(mode).alias('geometry'),
            pl.lit(population).alias('population'),pl.lit(source).alias('source_file'),
            pl.lit(records[source]['sha256']).alias('source_sha256'),pl.lit('verified_final_artifact').alias('source_status'))
    def read(name): return pl.read_parquet(root/name)
    generated = read('repertoire_metrics.parquet')
    if generated.height != 2400 or generated['chain_id'].n_unique()!=2400:
        raise ValueError('physical chain count/keys changed')
    fraction=CASES[case][0]
    memory=int(1200*fraction)
    expected_counts={(component,locus):count for component,n in [('naive',1200-memory),('memory',memory)]
                     for locus,count in [('IGH',n),('IGK',n//2),('IGL',n//2)]}
    actual_counts={(r['component'],r['locus']):r['len'] for r in generated.group_by('component','locus').len().to_dicts()}
    if actual_counts!=expected_counts or not generated['rna_pairs'].eq(100).all():
        raise ValueError('component/locus/RNA input counts changed')
    if generated.filter((pl.col('locus')=='IGH') & (pl.col('constant_nt_gene')!=
            pl.when(pl.col('component')=='naive').then(pl.lit('IGHM')).otherwise(pl.lit('IGHG1')))).height:
        raise ValueError('heavy constant assignments changed')
    scalars = generated.drop('variable_aa','protein_aa','generator_junction_nt','generator_junction_aa',
                              'annotated_junction_nt','annotated_junction_aa').with_columns(
        generated['generator_junction_aa'].str.len_chars().alias('junction_length_aa'),
        (generated['generator_junction_aa'].str.len_chars()-2).alias('generator_cdr3_length_aa'))
    output['generated_chains'] = tag(scalars,'repertoire_metrics.parquet','generated_chain_records')
    heavy = scalars.filter(pl.col('locus')=='IGH')
    if (heavy.height!=1200 or heavy['antibody_id'].n_unique()!=1200
            or heavy.filter(pl.col('component')=='memory').height!=int(1200*fraction)):
        raise ValueError('observed antibody/component denominator changed')
    antibodies=heavy.with_columns(pl.lit(1200).alias('antibody_denominator'),(pl.col('rna_pairs')/pl.col('rna_pairs').sum()).alias('rna_antibody_fraction'),
        (pl.col('serum_weight')/pl.col('serum_weight').sum()).alias('protein_antibody_fraction'))
    output['generated_antibodies']=tag(antibodies,'repertoire_metrics.parquet','generated_observed_antibody_records')
    candidates = read(mode+'_candidate_identity.parquet')
    candidate_scalars=candidates.drop('seq','nt_seq','cdr3_seq','region_ids').with_columns(
        candidates['seq'].str.len_chars().alias('observed_fragment_length_aa'),
        candidates['cdr3_seq'].str.len_chars().alias('observed_junction_length_aa'))
    output['observed_candidates']=tag(candidate_scalars,mode+'_candidate_identity.parquet','observed_heavy_candidate_records')
    with gzip.open(root/'assembly/arda_contigs.json.gz','rt') as stream: annotations=json.load(stream)
    annotation=pl.DataFrame([dict(sequence_id=a['sequence_id'],annotation_locus=a.get('locus'),
        v_call=a.get('v_call'),j_call=a.get('j_call'),c_call=a.get('c_call'),
        annotated_junction_length_aa=len(a['junction_aa']) if a.get('junction_aa') else None) for a in annotations],
        schema={'sequence_id':pl.String,'annotation_locus':pl.String,'v_call':pl.String,'j_call':pl.String,
                'c_call':pl.String,'annotated_junction_length_aa':pl.Int64})
    catalog=read('contig_catalog.parquet')
    if set(annotation['sequence_id'])!=set(catalog['sequence_id']): raise ValueError('annotation/catalog keys differ')
    observed=catalog.join(annotation,on='sequence_id',how='left',validate='1:1')
    output['observed_contigs']=tag(observed,'contig_catalog.parquet','observed_raw_contig_records')
    gene_tables=[gene_frequencies(generated,{'V':'v_call','J':'j_call','C':'constant_nt_gene'},'generated_chain_records','native_GenAIRR_OGRDB'),
        gene_frequencies(observed,{'V':'v_call','J':'j_call','C':'c_call'},'observed_raw_contig_records','native_ARDA_calls'),
        gene_frequencies(candidates.with_columns(pl.lit('IGH').alias('locus')),
                         {'V':'v_name','J':'j_name','C':'c_name'},'observed_heavy_candidate_records','native_ARDA_calls')]
    # Keep separate namespaces/populations and every missing/ambiguous call bin.
    for frame,source in zip(gene_tables,['repertoire_metrics.parquet','assembly/arda_contigs.json.gz',mode+'_candidate_identity.parquet']):
        frame=frame.with_columns(pl.lit(case).alias('case'),pl.lit(mode).alias('geometry'),pl.lit(source).alias('source_file'),
            pl.lit(records[source]['sha256']).alias('source_sha256'),pl.lit('verified_final_artifact').alias('source_status'))
        output.setdefault('gene_usage',[]).append(frame)
    output['gene_usage']=pl.concat(output['gene_usage'],how='vertical_relaxed')
    for name in EXISTING:
        frame=read(mode+'_'+name+'.parquet')
        if name=='identity_population': frame=frame.drop([c for c in ('variable_aa','protein_aa','generator_junction_nt','generator_junction_aa','annotated_junction_nt','annotated_junction_aa') if c in frame.columns])
        output[name]=tag(frame,mode+'_'+name+'.parquet',
            'physical_heavy_records' if name=='identity_population' else
            'producer_defined_heavy_identity_readouts' if name=='identity_metrics' else
            'producer_defined_pooled_coverage_or_recovery')
    output['component_counts']=tag(read('component_counts.parquet'),'component_counts.parquet','generated_chain_records')
    calls=[];native=[];species=[]
    for enzyme in ('Trypsin','Chymotrypsin'):
        name=mode+'_'+enzyme+'_psm_correctness.parquet'
        calls.append(read(name))
        native.append(tag(read(mode+'_'+enzyme+'_readout_metrics.parquet').with_columns(pl.lit(enzyme).alias('enzyme')),
            mode+'_'+enzyme+'_readout_metrics.parquet','native_PSM_assignments_and_events'))
    for enzyme,frame in [('Trypsin',calls[0]),('Chymotrypsin',calls[1]),('pooled',pl.concat(calls,how='vertical_relaxed'))]:
        n=eligible_species(frame).height
        species.append(dict(case=case,geometry=mode,enzyme=enzyme,quantity='canonical_IL_called_species_ge6aa',
            count=n,numerator=n,denominator=None,n_observations=frame.height,unit='species',population='called_PSM_assignment_records',
            source_status='verified_final_artifacts',
            source_file_1=mode+'_'+('Trypsin' if enzyme=='pooled' else enzyme)+'_psm_correctness.parquet',
            source_sha256_1=records[mode+'_'+('Trypsin' if enzyme=='pooled' else enzyme)+'_psm_correctness.parquet']['sha256'],
            source_file_2=mode+'_Chymotrypsin_psm_correctness.parquet' if enzyme=='pooled' else None,
            source_sha256_2=records[mode+'_Chymotrypsin_psm_correctness.parquet']['sha256'] if enzyme=='pooled' else None))
    output['native_readouts']=pl.concat(native,how='vertical_relaxed')
    pooled=pl.concat(calls,how='vertical_relaxed')
    evaluated=pooled.filter(pl.col('correct').is_not_null())
    correct=int(evaluated['correct'].sum() or 0)
    output['pooled_PSM_accuracy']=pl.DataFrame([dict(case=case,geometry=mode,enzyme='pooled',
        correct_assignments=correct,evaluated_assignments=evaluated.height,
        all_assignments=pooled.height,unscorable_assignments=pooled.height-evaluated.height,
        numerator=correct,denominator=evaluated.height,n_observations=pooled.height,unit='assignment_fraction',
        fraction=correct/evaluated.height if evaluated.height else None,
        source_file_1=mode+'_Trypsin_psm_correctness.parquet',
        source_sha256_1=records[mode+'_Trypsin_psm_correctness.parquet']['sha256'],
        source_file_2=mode+'_Chymotrypsin_psm_correctness.parquet',
        source_sha256_2=records[mode+'_Chymotrypsin_psm_correctness.parquet']['sha256'],
        source_status='verified_final_artifacts')])
    output['peptide_species']=pl.DataFrame(species)
    return output


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--cpus',type=int)
    args=parser.parse_args(argv)
    if args.output.exists(): parser.error('output must be fresh')
    resolve_cpu_budget(args.cpus)
    import polars as pl
    inputs=inventory(args.input_root.resolve())
    validated=[(case,mode,root,*validate_bundle(root,case,mode)) for case,mode,root in inputs]
    tables={}; provenance=[]
    for case,mode,root,records,report,monitor in validated:
        for name,frame in summarize_tables(root,case,mode,records).items(): tables.setdefault(name,[]).append(frame)
        provenance.append(dict(case=case,geometry=mode,bundle_manifest_sha256=digest(root/'retained_table_manifest.json'),
            report_sha256=records['metric_integration.json']['sha256'],producer_software_sha256=report['software_sha256'],
            monitoring_source_sha256=monitor,
            sources={name:r['sha256'] for name,r in records.items()}))
    # Validate and form all tables before publishing any output.
    final={}
    for name,parts in sorted(tables.items()):
        frame=pl.concat(parts,how='vertical_relaxed')
        final[name]=frame.sort(frame.columns,nulls_last=True)
    args.output.mkdir(parents=True)
    for name,frame in final.items():
        with (args.output/(name+'.tsv.gz')).open('wb') as raw:
            with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as stream:
                frame.write_csv(stream,separator='\t',null_value='null')
    (args.output/'provenance.json').write_text(json.dumps(dict(cases=provenance,
        output_sha256={p.name:digest(p) for p in sorted(args.output.glob('*.tsv.gz'))},
        units={'lengths':'amino-acid residues, unless named nt', 'branch_or_ancestor_substitutions':'nucleotide endpoint differences, not summed independent mutation events',
               'RNA':'read pairs', 'serum_weight':'declared native protein input units',
               'fractions':'dimensionless; inherited chain-pool fractions differ from explicitly renormalized antibody fractions',
               'counts':'exact records/classes/events as individually named; never donors or independent replicates'},
        summary_script_sha256=digest(Path(__file__)),coverage_scope='pooled across proteases; existing three evidence levels; overlapping views never summed',
        inference='no tests, controls, thresholds, scoring, fitting or automatic selection',
        distributions='unbinned scalar chain/antibody records; no inferred clone uniqueness or biological replicates',
        gene_namespaces='native generation and native observed annotation bins remain distinct; canonical V fields preserved'),indent=2)+'\n')


if __name__=='__main__': main()
