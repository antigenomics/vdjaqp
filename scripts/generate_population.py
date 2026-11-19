"""Generate a conditional paired fixed-heavy-constant population using native GenAIRR only.

Resource qualification/mechanistic substrate; not a fitted experimental mock.
"""
import argparse
import hashlib
from importlib.resources import files
import json
import math
from pathlib import Path
import shutil
import sys
from importlib.metadata import version
import numpy as np
from Bio.Seq import translate
from GenAIRR import Experiment

sys.path.insert(0, str(Path(__file__).resolve().parent))
from smoke_generate import load_constants, reference_hashes
from validate_constant_candidates import load_validated_constants, MANIFEST as CONSTANT_MANIFEST

from genairr_gene_map import load_reference_gene_mapping, canonical_v_fields, write_reference_gene_mapping


HEAVY_CONSTANT_GENES = ('IGHG1', 'IGHM', 'IGHG3', 'IGHA1', 'IGHE')

def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('duplicate model key: '+key)
        result[key] = value
    return result

def load_recombination_parameters(path):
    """Validate a complete heavy-only native model; never infer missing weights."""
    if path is None:
        return {}, None
    raw = Path(path).read_bytes()
    model = json.loads(raw, object_pairs_hook=unique_object)
    if not isinstance(model, dict) or set(model) != {'schema_version', 'genairr_version', 'reference_sha256', 'parameters'}:
        raise ValueError('invalid recombination model fields')
    identity = json.loads((Path(__file__).resolve().parents[1]/'envs/generation/identity.json').read_text())
    reference = files('GenAIRR.data')/'builtin_dataconfigs/HUMAN_IGH_OGRDB.pkl'
    digest = hashlib.sha256(reference.read_bytes()).hexdigest()
    if (type(model['schema_version']) is not int or model['schema_version'] != 1 or
            model['genairr_version'] != version('GenAIRR') or
            model['genairr_version'] != identity['versions']['GenAIRR'] or
            model['reference_sha256'] != digest or
            digest != identity['genairr_reference_sha256']['HUMAN_IGH_OGRDB.pkl']):
        raise ValueError('recombination model version/reference mismatch')
    from GenAIRR.data import HUMAN_IGH_OGRDB as config
    params = model['parameters']
    if not isinstance(params, dict) or set(params) != {'v_allele_weights', 'j_allele_weights', 'np1_lengths', 'np2_lengths'}:
        raise ValueError('only complete V/J weights and NP1/NP2 lengths are supported')
    def positive(value):
        return type(value) in (int, float) and 0 < value <= sys.float_info.max and math.isfinite(value)
    for key, alleles in [('v_allele_weights', config.v_alleles), ('j_allele_weights', config.j_alleles)]:
        names = {allele.name for group in alleles.values() for allele in group}
        weights = params[key]
        if not isinstance(weights, dict) or set(weights) != names or not all(positive(w) for w in weights.values()):
            raise ValueError('require every native allele with finite positive weight: '+key)
    for key in ['np1_lengths', 'np2_lengths']:
        pairs = params[key]
        if (not isinstance(pairs, list) or len(pairs) != 41 or
                any(not isinstance(p, list) or len(p) != 2 or type(p[0]) is not int or not positive(p[1]) for p in pairs) or
                {p[0] for p in pairs} != set(range(41))):
            raise ValueError('require original NP length support 0..40 with finite positive weights: '+key)
    return params, dict(schema_version=1, sha256=hashlib.sha256(raw).hexdigest(), model=model)


def draw_seed(seed, stage, locus, family, node, attempt):
    """Injective mixed-radix coordinates; root batch reserves [seed, seed+families)."""
    bounds = ((seed, 10**14), (stage, 3), (locus, 3),
              (family, 10**10), (node, 3), (attempt, 100))
    if any(type(value) is not int or not 0 <= value < limit for value, limit in bounds) or stage == 0:
        raise ValueError('seed coordinates outside declared mixed-radix bounds')
    return seed + stage * 10**15 + locus * 10**14 + family * 1000 + node * 100 + attempt


def compatibility_reason(record, constant):
    nt, aa = record['sequence'].upper(), record['sequence_aa']
    if not record['productive']:
        return 'nonproductive'
    if str(translate(nt[:len(nt)//3*3])) != aa:
        return 'variable_translation_mismatch'
    full = nt + constant['nt']
    if len(full) % 3:
        return 'constant_frame_mismatch'
    if str(translate(full)) != aa + constant['aa'] + '*':
        return 'constant_translation_mismatch'
    return None


def distance(nt, ancestor):
    if len(nt) != len(ancestor):
        raise ValueError('native continuation changed sequence length')
    if not nt:
        raise ValueError('native sequence is empty')
    n = int(np.count_nonzero(np.frombuffer(nt.encode('ascii'), dtype=np.uint8) !=
                             np.frombuffer(ancestor.encode('ascii'), dtype=np.uint8)))
    return {'substitutions': n, 'denominator_nt': len(nt), 'fraction': n/len(nt)}


def load_population_constants(reference_dir, heavy_constant_gene='IGHG1', additional_reference_dir=None):
    """Select exact validated source suffixes; retain the historical default loader."""
    if heavy_constant_gene not in HEAVY_CONSTANT_GENES:
        raise ValueError('unknown heavy constant gene')
    if heavy_constant_gene == 'IGHG1' and additional_reference_dir is not None:
        raise ValueError('additional-reference-dir is only valid for a nondefault heavy constant')
    if heavy_constant_gene != 'IGHG1' and additional_reference_dir is None:
        raise ValueError('nondefault heavy constant requires additional-reference-dir')
    reference_dir = Path(reference_dir)
    expected = {}
    manifest_dir = Path(__file__).resolve().parents[1]/'envs/generation'
    for name in ('references.json', 'lambda_references.json'):
        expected.update(json.loads((manifest_dir/name).read_text())['references'])
    hashes = reference_hashes(reference_dir, all_loci=True)
    for name, metadata in expected.items():
        if hashes[str(reference_dir/name)] != metadata['sha256']:
            raise ValueError(f'raw reference hash mismatch: {name}')
    constants = load_constants(reference_dir, all_loci=True)
    additional = None
    if additional_reference_dir is not None:
        # Validate every pinned record, including unused candidates, before any write/draw.
        selected = load_validated_constants(additional_reference_dir)
        constants['IGH'] = selected[heavy_constant_gene]
        from validate_constant_candidates import validated_source_manifest
        additional = validated_source_manifest(additional_reference_dir, CONSTANT_MANIFEST)
    return constants, expected, additional


def write_constant_references(path, constants):
    """Stage declared public constant proteins as search model inputs, never chain truth."""
    Path(path).write_text(''.join(
        f">constant_{constants[locus]['protein_accession']}\n{constants[locus]['aa']}\n"
        for locus in constants))


def mixture_layout(total_antibodies, memory_fraction=None, protein_concentration='equal', memory_records=None):
    """Exact observed counts and component/light-balanced head-family membership."""
    if type(total_antibodies) is not int or total_antibodies <= 0:
        raise ValueError('total-antibodies must be a positive integer')
    if memory_records is not None:
        if memory_fraction is not None or type(memory_records) is not int or not 0 <= memory_records <= total_antibodies:
            raise ValueError('memory-records must be an integer within total, alternative to memory-fraction')
        memory_count = memory_records
    else:
        if memory_fraction not in (.25, .5, .75):
            raise ValueError('memory-fraction must be 0.25, 0.5 or 0.75')
        memory_count = total_antibodies * memory_fraction
    if memory_count != int(memory_count) or int(memory_count) % 6:
        raise ValueError('memory antibody count must be divisible by six for balanced three-node trees')
    memory_count = int(memory_count)
    naive_count = total_antibodies - memory_count
    if naive_count % 2:
        raise ValueError('naive antibody count must be even for balanced light loci')
    memory_families = memory_count // 3
    if naive_count + memory_families >= 10**10:
        raise ValueError('mixture exceeds family seed bounds')
    heads = set()
    if protein_concentration == 'head80':
        if naive_count % 20 or memory_families % 20:
            raise ValueError('head80 requires each component family count divisible by 20 for exact light-balanced ten percent')
        for start, count in ((0, naive_count), (naive_count, memory_families)):
            for light in (0, 1):
                heads.update(range(start + light, start + count // 10, 2))
    return dict(naive_antibodies=naive_count, memory_antibodies=memory_count,
                memory_families=memory_families, families=naive_count + memory_families,
                head_families=heads)


def load_heavy_allocation(path, mixture):
    """Require one declared constant per observed antibody, independently of SHM."""
    raw = Path(path).read_bytes()
    model = json.loads(raw, object_pairs_hook=unique_object)
    if (not isinstance(model, dict) or set(model) != {'schema_version', 'allocation'} or
            type(model['schema_version']) is not int or model['schema_version'] != 1):
        raise ValueError('invalid heavy constant allocation schema')
    expected = {f'family{family}_node{node}' for family in range(mixture['families'])
                for node in range(1 if family < mixture['naive_antibodies'] else 3)}
    allocation = model['allocation']
    if (not isinstance(allocation, dict) or set(allocation) != expected or
            any(type(gene) is not str or gene not in HEAVY_CONSTANT_GENES for gene in allocation.values())):
        raise ValueError('heavy allocation must name exactly every observed antibody with a supported gene')
    return allocation, dict(sha256=hashlib.sha256(raw).hexdigest(), model=model)


def generate(output, reference_dir, families=None, pairs_per_chain=100, mutation_rate=.01,
             protein_concentration='equal', seed=20261110, heavy_constant_gene='IGHG1',
             additional_reference_dir=None, total_antibodies=None, memory_fraction=None,
             transcript_reference_dir=None, transcript_threads=16, recombination_parameters=None,
             memory_records=None, heavy_constant_allocation=None):
    mixture = None
    if total_antibodies is not None or memory_fraction is not None or memory_records is not None:
        if families is not None:
            raise ValueError('choose families or total-antibodies/memory-fraction, not both')
        if heavy_constant_gene != 'IGHG1' or additional_reference_dir is None:
            raise ValueError('mixture requires memory IGHG1 and explicit additional-reference-dir for naive IGHM')
        mixture = mixture_layout(total_antibodies, memory_fraction, protein_concentration, memory_records)
        families = mixture['families']
    if type(families) is not int or not 0 < families < 10**10 or families % 2:
        raise ValueError('families must be positive, even and below 10000000000')
    if type(pairs_per_chain) is not int or pairs_per_chain <= 0:
        raise ValueError('pairs-per-chain must be a positive integer')
    if not math.isfinite(mutation_rate) or not 0 <= mutation_rate <= 1:
        raise ValueError('mutation-rate must be finite in [0,1]')
    if type(seed) is not int or not 0 <= seed < 10**14:
        raise ValueError('seed must be an integer in [0,100000000000000)')
    if protein_concentration not in ('equal', 'head80'):
        raise ValueError('unknown protein concentration')
    if type(transcript_threads) is not int or transcript_threads < 1:
        raise ValueError('transcript-threads must be positive')
    output, reference_dir = Path(output), Path(reference_dir)
    if output.exists():
        raise FileExistsError('output must be a fresh absent directory')
    allocation = allocation_provenance = heavy_constants = None
    if heavy_constant_allocation is not None:
        if mixture is None or additional_reference_dir is None:
            raise ValueError('heavy allocation requires explicit mixture and additional-reference-dir')
        allocation, allocation_provenance = load_heavy_allocation(heavy_constant_allocation, mixture)
    heavy_recombination, recombination_provenance = load_recombination_parameters(recombination_parameters)
    constants, expected, additional = load_population_constants(
        reference_dir, 'IGHM' if mixture else heavy_constant_gene, additional_reference_dir)
    if mixture:
        naive_constant = constants['IGH']
        constants['IGH'] = load_constants(reference_dir, all_loci=True)['IGH']
        constants['IGH_naive'] = naive_constant
    if allocation is not None:
        heavy_constants = load_validated_constants(additional_reference_dir)
        heavy_constants['IGHG1'] = constants['IGH']
    manifest_dir = Path(__file__).resolve().parents[1]/'envs/generation'
    gene_mapping = load_reference_gene_mapping()
    output.mkdir(parents=True)
    (output/'references').mkdir()
    write_reference_gene_mapping(output/'reference_gene_mapping.tsv', gene_mapping)
    for name in expected:
        shutil.copyfile(reference_dir/name, output/'references'/name)
    if additional is not None:
        target = output/'references/additional_constants'
        target.mkdir()
        for name in additional['downloaded_records']:
            shutil.copyfile(Path(additional_reference_dir)/name, target/name)
        (target/'manifest.json').write_text(json.dumps(additional, indent=2)+'\n')
    write_constant_references(output/'constant_reference_proteins.fasta',
        dict(heavy_constants, IGK=constants['IGK'], IGL=constants['IGL']) if allocation is not None else constants)
    transcript_contexts = None
    if transcript_reference_dir is not None:
        import csv
        import subprocess
        from build_transcript_context import prepare
        source_model = output/'references/transcript_model_sources.json'
        prepare(transcript_reference_dir, source_model)
        repo = Path(__file__).resolve().parents[1]
        subprocess.run([str(repo/'.venv/bin/python'), str(repo/'scripts/build_transcript_context.py'),
                        '--align', str(source_model), '--output', str(output/'references/transcript_model'),
                        '--threads', str(transcript_threads)], check=True, timeout=300)
        with (output/'references/transcript_model/contexts.tsv').open() as stream:
            transcript_contexts = {r['native_allele']: r for r in csv.DictReader(stream, delimiter='\t')}
        # Generation has Bio/GenAIRR, not the scoring package's native toolchain.
        import importlib.util
        spec = importlib.util.spec_from_file_location('transcript_context', repo/'vdjaqp/transcript_context.py')
        context_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(context_module)
        leader_rna_template = context_module.leader_rna_template
    rows, attempts, latent_ancestors = [], [], []
    head = len(mixture['head_families']) if mixture else math.ceil(.1*families)
    for li, (locus, cfg, model, gene) in enumerate([
            ('IGH', 'human_igh', 'hh_s5f', heavy_constant_gene),
            ('IGK', 'human_igk', 'hkl_s5f', 'IGKC'),
            ('IGL', 'human_igl', 'hkl_s5f', 'IGLC2')]):
        family_ids = list(range(families)) if locus == 'IGH' else list(range(li-1, families, 2))
        exp = Experiment.on(cfg)
        if locus == 'IGL':
            exp = exp.restrict_alleles(j=['IGLJ2/3-HDWH*01','IGLJ2/3-JWS7*02'])
        base = exp.recombine(**(heavy_recombination if locus == 'IGH' else {})).productive_only().compile()
        mutator = Experiment.on(cfg).mutate(model='s5f', rate=mutation_rate, s5f_model=model).productive_only().compile()
        root_seed = seed + li * 10**14
        roots = base.run(n=len(family_ids), seed=root_seed, strict=True)
        for index, family in enumerate(family_ids):
            component = ('naive' if family < mixture['naive_antibodies'] else 'memory') if mixture else None
            constant = constants['IGH_naive'] if locus == 'IGH' and component == 'naive' else constants[locus]
            selected_gene = 'IGHM' if locus == 'IGH' and component == 'naive' else gene
            ancestor_id = f'family{family}_ancestor_{locus}' if component == 'memory' else None
            accepted = []
            nodes = (-1, 0, 1, 2) if component == 'memory' else range(1 if component == 'naive' else 3)
            for node in nodes:
                recombining = node == -1 or (node == 0 and component != 'memory')
                if allocation is not None and locus == 'IGH':
                    selected_gene = allocation[f'family{family}_node{max(node, 0)}']
                    constant = heavy_constants[selected_gene]
                for attempt in range(100):
                    if recombining and attempt == 0:
                        native, used_seed = roots[index], root_seed + index
                    else:
                        used_seed = draw_seed(seed, 1 if recombining else 2, li, family, max(node, 0), attempt)
                        native = (base.run(n=1, seed=used_seed, strict=True)[0] if recombining else
                                  mutator.simulator.run_from(accepted[-1].final_simulation(), used_seed, strict=True))
                    record = native.airr_record(base.refdata)
                    reason = compatibility_reason(record, constant)
                    counters = {k: v for k, v in record.items() if k.startswith('n_') or k == 'mutation_rate'}
                    attempts.append(dict(family_id=f'family{family}', locus=locus, node=node,
                                         attempt=attempt, seed=used_seed, rejection_reason=reason,
                                         native_counters=counters))
                    if reason is None:
                        break
                else:
                    (output/'failed_attempts.json').write_text(json.dumps(attempts, indent=2)+'\n')
                    raise RuntimeError(f'constant compatibility exhausted 100 native draws: {family}/{locus}/{node}')
                nt, aa = record['sequence'].upper(), record['sequence_aa']
                if node == -1:
                    accepted.append(native)
                    latent_ancestors.append(dict(chain_id=ancestor_id, family_id=f'family{family}',
                        locus=locus, component='memory', observed=False, variable_nt=nt, variable_aa=aa,
                        accepted_seed=used_seed, airr_truth=record, constant_nt_gene=selected_gene))
                    continue
                parent_nt = accepted[-1].airr_record(base.refdata)['sequence'].upper() if accepted else nt
                root_nt = accepted[0].airr_record(base.refdata)['sequence'].upper() if accepted else nt
                accepted.append(native)
                intensity = 1000.0
                if protein_concentration == 'head80':
                    if mixture:
                        intensity *= 8 if family in mixture['head_families'] else 2/9
                    else:
                        intensity *= .8*families/head if family < head else .2*families/(families-head)
                coding = nt[:len(nt)//3*3]
                rows.append(dict(chain_id=f'family{family}_node{node}_{locus}',
                    antibody_id=f'family{family}_node{node}', family_id=f'family{family}',
                    parent_id=ancestor_id if node == 0 else f'family{family}_node{node-1}_{locus}',
                    locus=locus, variable_nt=nt, coding_variable_nt=coding,
                    untranslated_terminal_nt=nt[len(coding):], variable_aa=aa,
                    protein_aa=aa+constant['aa'], constant_accession=constant['protein_accession'],
                    constant_cds=constant['cds_accession'], constant_nt_gene=selected_gene,
                    amplicon_nt=nt+constant['nt'][:29], rnaseq_nt=nt+constant['nt'],
                    rna_pairs=pairs_per_chain, secretion_factor=intensity/pairs_per_chain,
                    protein_intensity=intensity, branch_distance=distance(nt,parent_nt),
                    root_distance=distance(nt,root_nt), accepted_seed=used_seed, airr_truth=record,
                    **canonical_v_fields(locus, record['v_call'], gene_mapping)))
                if transcript_contexts is not None:
                    context = transcript_contexts[record['v_call']]
                    template = leader_rna_template(nt, constant['nt'], aa+constant['aa'], context,
                        v_trim_5=record.get('v_trim_5'), translate_nt=lambda seq: str(translate(seq)))
                    rows[-1].update(rnaseq_nt=template['rnaseq_nt'], rna_model='leader_imputed_v1',
                        rna_mature_nt_start=template['mature_nt_start'],
                        rna_leader_nt_start=template['leader_nt_start'], rna_prefix_nt=context['prefix_nt'],
                        rna_leader_nt=context['leader_nt'], rna_context_method=context['method'],
                        rna_context_confidence=context['confidence'],
                        rna_context_source_id=context['selected_source_id'])
                if mixture:
                    rows[-1].update(component=(('unmutated_singleton' if component == 'naive' else 'related_shm_family')
                                              if allocation is not None else component), ancestor_id=ancestor_id,
                                    is_head=family in mixture['head_families'])
    rows.sort(key=lambda r: (int(r['family_id'][6:]), int(r['antibody_id'].split('_node')[1]), r['locus']))
    for name, field in [('proteins.fasta','protein_aa'),('amplicon_templates.fasta','amplicon_nt'),('rnaseq_templates.fasta','rnaseq_nt')]:
        (output/name).write_text(''.join(f'>{r["chain_id"]}'+(f' [# intensity={r["protein_intensity"]} #]' if name=='proteins.fasta' else '')+f'\n{r[field]}\n' for r in rows))
    (output/'read_counts.tsv').write_text(''.join(f'{r["chain_id"]}\t{2*r["rna_pairs"]}\n' for r in rows))
    truth = dict(purpose=f"conditional native {'IgG1' if heavy_constant_gene == 'IGHG1' else heavy_constant_gene} resource/mechanistic population; not fitted reconstruction",
        seed=seed, chains=rows, reference_gene_mapping=gene_mapping, reference_sha256={name: m['sha256'] for name,m in expected.items()},
        constant_reference_proteins='constant_reference_proteins.fasta',
        constant_references=constants, family_topology='root-child-grandchild',
        lambda_scope='native ambiguous J2/3 pool restricted to fixed IGLC2; balanced kappa/lambda assumption',
        parameters=dict(families=families,pairs_per_chain=pairs_per_chain,mutation_rate=mutation_rate,
                        protein_concentration=protein_concentration, head_family_count=head,
                        head_order='ascending family ID', max_attempts=100, heavy_constant_gene=heavy_constant_gene),
        conditioning='whole native draws rejected for incompatibility with fixed selected constant protein; no repair or J masking',
        attempts=attempts, rejected_draws=sum(a['rejection_reason'] is not None for a in attempts),
        runtime='GenAIRR native run_batch serial map; no parallel wrapper',
        versions={name: version(name) for name in ('GenAIRR','biopython')},
        engine_identity=json.loads((manifest_dir/'identity.json').read_text()),
        seed_policy='root batches seed+locus*10^14+i; retries/continuations injective stage/locus/family/node/attempt mixed radix',
        abundance_scope='externally prescribed equal RNA pairs; arbitrary molar protein units; RNA pairs times secretion factor equals intensity; paired H/L equal, common 2H2L factor uncalibrated',
        constant_nt_scope='full annotated spliced C CDS plus native variable; terminal stop included; no leader or UTR')
    if recombination_provenance is not None:
        truth['heavy_recombination_model'] = recombination_provenance
        truth['producer_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if mixture:
        truth['purpose'] = 'conditional native naive IgM / memory IgG1 mixture; not fitted reconstruction'
        truth['family_topology'] = 'independent observed naive singleton; latent memory ancestor to three observed descendants'
        truth['latent_ancestors'] = latent_ancestors
        truth['component_counts'] = {k: mixture[k] for k in
            ('naive_antibodies', 'memory_antibodies', 'memory_families')}
        truth['component_counts'].update(observed_antibodies=total_antibodies,
                                         observed_chains=len(rows), latent_ancestor_chains=len(latent_ancestors))
        truth['parameters'].update(total_antibodies=total_antibodies, memory_fraction=memory_fraction,
                                  naive_constant_gene='IGHM', memory_constant_gene='IGHG1',
                                  head_order='ascending family ID within component and light locus')
        truth['abundance_scope'] += '; head80 selects exact ten percent families within each component/light stratum; 80 percent component protein mass; RNA pairs remain equal per observed chain'
    if memory_records is not None:
        truth['parameters']['memory_records'] = memory_records
    if allocation is not None:
        truth['purpose'] = 'conditional native constant-allocated reconstruction substrate; construction components are not phenotypes'
        truth['heavy_constant_allocation'] = allocation_provenance
        truth['heavy_constant_references'] = heavy_constants
        truth['producer_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        truth['parameters'].pop('naive_constant_gene', None)
        truth['parameters'].pop('memory_constant_gene', None)
        truth['family_topology'] = 'independent unmutated observed singleton; latent ancestor to three observed SHM descendants'
        truth['conditioning'] += '; allocated constant selected per observed node; latent ancestor conditioned on first descendant constant only, never assayed'
    if additional is not None:
        truth['constant_nt_scope'] = ('annotated spliced C CDS suffix plus unmodified native variable; '
                                      'overlap_one_identical_G removes only duplicate reference G; '
                                      'terminal stop included; no leader or UTR')
        truth['additional_constant_sources'] = additional
        truth['reference_sha256'].update({
            'additional_constants/'+name: source['sha256']
            for name, source in additional['downloaded_records'].items()})
    if transcript_contexts is not None:
        truth['rna_model'] = 'leader_imputed_v1'
        truth['transcript_model'] = json.loads((output/'references/transcript_model/audit.json').read_text())
        truth['constant_nt_scope'] += '; RNA-only validated/donor-imputed spliced leader prefix; no fabricated UTR/poly(A); mature protein and amplicon unchanged'
    (output/'truth.json').write_text(json.dumps(truth, indent=2)+'\n')
    return rows


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True)
    p.add_argument('--reference-dir', required=True)
    p.add_argument('--families', type=int)
    p.add_argument('--total-antibodies', type=int)
    p.add_argument('--memory-records', type=int, help='exact related-family record count, alternative to memory-fraction')
    p.add_argument('--heavy-constant-allocation', help='schema-1 exact antibody-ID to heavy-gene JSON allocation')
    p.add_argument('--memory-fraction', type=float, choices=[.25, .5, .75])
    p.add_argument('--pairs-per-chain', type=int, default=100)
    p.add_argument('--mutation-rate', type=float, default=.01)
    p.add_argument('--protein-concentration', choices=['equal','head80'], default='equal')
    p.add_argument('--heavy-constant-gene', choices=HEAVY_CONSTANT_GENES, default='IGHG1')
    p.add_argument('--additional-reference-dir')
    p.add_argument('--seed', type=int, default=20261110)
    p.add_argument('--transcript-reference-dir', help='enable separately versioned leader-imputed RNA model from pinned raw exports')
    p.add_argument('--recombination-parameters', help='versioned heavy-only native V/J and NP model JSON')
    p.add_argument('--transcript-threads', type=int, default=16)
    args = p.parse_args()
    print('Generated', len(generate(**vars(args))), 'chains')
