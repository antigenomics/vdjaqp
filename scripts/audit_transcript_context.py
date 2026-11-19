"""Audit pinned public transcript-context parameters, without drawing a repertoire.

Rebuild requested coverage tables from raw references; never feed audit outputs
back as a generator cache. Trusted native pickle bytes are checked before loading.
"""
import argparse
from collections import Counter
import csv
from decimal import Decimal
import hashlib
from importlib.resources import files
import json
from pathlib import Path

from Bio.Seq import translate
from genairr_gene_map import exact_v_map

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO/'envs/generation/transcript_context_sources.json'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def normalized(sequence):
    return sequence.upper().replace('.', '').replace('-', '')


def extract_context(allele):
    """One-based closed source coordinates; remove the L1/L2 intron only."""
    raw = allele['sequence'].upper()
    coding = normalized(allele['coding_sequence'])
    start, end = allele.get('gene_start'), allele.get('gene_end')
    if (type(start) is not int or type(end) is not int
            or not 1 <= start <= end <= len(raw)
            or raw[start-1:end] != coding):
        return dict(status='coding_source_coordinate_mismatch')
    keys = ('leader_1_start','leader_1_end','leader_2_start','leader_2_end')
    coords = [allele.get(k) for k in keys]
    if any(value is None for value in coords):
        return dict(status='missing_leader_coordinates')
    l1, r1, l2, r2 = coords
    if (any(type(v) is not int for v in coords)
            or not 1 <= l1 <= r1 < l2 <= r2 < start <= end <= len(raw)
            or r2 + 1 != start):
        return dict(status='invalid_leader_coordinates')
    leader = raw[l1-1:r1] + raw[l2-1:r2]
    if (set(leader+coding)-set('ACGT') or not leader.startswith('ATG')
            or len(leader) % 3 or '*' in str(translate(leader))):
        return dict(status='invalid_leader_frame_start_or_stop',leader_length_nt=len(leader),
                    leader_start_atg=leader.startswith('ATG'),leader_frame_zero=len(leader)%3==0,
                    leader_stop_free='*' not in str(translate(leader[:len(leader)//3*3])))
    # Frame-zero core projection must survive the spliced leader junction.
    if str(translate(leader+coding[:len(coding)//3*3])) != str(translate(leader)) + str(translate(coding[:len(coding)//3*3])):
        return dict(status='leader_core_translation_mismatch')
    u1, u2 = allele.get('utr_5_prime_start'), allele.get('utr_5_prime_end')
    utr = ''
    if u1 is not None or u2 is not None:
        if (type(u1) is not int or type(u2) is not int
                or not 1 <= u1 <= u2 < l1 or u2+1 != l1):
            return dict(status='invalid_utr_coordinates')
        utr = raw[u1-1:u2]
        if set(utr)-set('ACGT'):
            return dict(status='invalid_utr_sequence')
    prefix = utr+leader
    return dict(status='valid', leader_nt=leader, leader_length_nt=len(leader),
                utr_nt=utr, utr_length_nt=len(utr), prefix_nt=prefix,
                prefix_sha256=hashlib.sha256(prefix.encode()).hexdigest())


def self_check():
    """Fixed nucleotide bookkeeping check; no stochastic or biological generation."""
    allele=dict(sequence='ATGGTCAGGCTCAG',coding_sequence='CAG',gene_start=12,gene_end=14,
                leader_1_start=1,leader_1_end=3,leader_2_start=9,leader_2_end=11)
    context=extract_context(allele)
    assert context['status']=='valid' and context['leader_nt']=='ATGGCT'
    assert context['prefix_nt']=='ATGGCT' and context['utr_length_nt']==0
    assert extract_context(dict(allele,gene_end=13))['status']=='coding_source_coordinate_mismatch'
    assert extract_context(dict(allele,leader_1_start=None))['status']=='missing_leader_coordinates'
    assert extract_context(dict(allele,sequence='TAAGTCAGGCTCAG'))['status']=='invalid_leader_frame_start_or_stop'
    assert extract_context(dict(allele,utr_5_prime_start=1))['status']=='invalid_utr_coordinates'


def load_verified_sources(raw_dir):
    """Validate pinned public and native reference bytes before trusted loading."""
    raw_dir = Path(raw_dir)
    specification = json.loads(MANIFEST.read_text())
    identity = json.loads((REPO/'envs/generation/identity.json').read_text())
    sources = {}
    # Verify ALL raw bytes/metadata and ALL pickle bytes before any pickle load.
    for locus, source in specification['sources'].items():
        payload = (raw_dir/source['filename']).read_bytes()
        require(hashlib.sha256(payload).hexdigest()==source['sha256'], f'{locus}: raw SHA256 mismatch')
        data = json.loads(payload, parse_float=Decimal)
        require(len(data['GermlineSet'])==1, f'{locus}: expected one germline set')
        germline = data['GermlineSet'][0]
        require(germline['locus']==locus and germline['species']['id']=='NCBITAXON:9606'
                and Decimal(str(germline['release_version']))==Decimal(str(source['release_version']))
                and germline['germline_set_ref']==source['germline_set_ref'], f'{locus}: source metadata mismatch')
        sources[locus] = [a for a in germline['allele_descriptions'] if a['sequence_type']=='V']
        require(all(a['locus']==locus and a['species']['id']=='NCBITAXON:9606'
                    for a in sources[locus]), f'{locus}: allele locus/species mismatch')
    base = files('GenAIRR.data').joinpath('builtin_dataconfigs')
    pickle_bytes = {}
    for locus in sources:
        name = f'HUMAN_{locus}_OGRDB.pkl'
        payload = base.joinpath(name).read_bytes()
        require(hashlib.sha256(payload).hexdigest()==identity['genairr_reference_sha256'][name],
                f'{locus}: installed native reference SHA256 mismatch')
        pickle_bytes[locus] = payload
    return specification, identity, sources, pickle_bytes


def audit(raw_dir, output):
    output = Path(output)
    require(not output.exists(), 'audit output must be fresh and absent')
    specification, identity, sources, pickle_bytes = load_verified_sources(raw_dir)
    import pickle
    coverage, matches, raw_rows, summaries = [], [], [], {}
    for locus, alleles in sources.items():
        contexts = {a['label']:extract_context(a) for a in alleles}
        require(len(contexts)==len(alleles), f'{locus}: duplicate source allele label')
        for allele in alleles:
            raw_rows.append(dict(locus=locus,source_allele=allele['label'],
                **{k:allele.get(k) for k in ('gene_start','gene_end','leader_1_start','leader_1_end',
                    'leader_2_start','leader_2_end','utr_5_prime_start','utr_5_prime_end')},
                **contexts[allele['label']]))
        config = pickle.loads(pickle_bytes[locus])
        legacy = [(a.name,a.gene,normalized(a.ungapped_seq)) for group in config.v_alleles.values() for a in group]
        canonical = [(a['label'],a['label'].split('*')[0],normalized(a['coding_sequence'])) for a in alleles]
        joined = exact_v_map(legacy, canonical)
        for name, row in joined.items():
            valid = {}; unresolved = []
            for match in row['canonical_matches']:
                label = match['allele']; context = contexts[label]
                matches.append(dict(locus=locus,native_allele=name,source_allele=label,**context))
                if context['status']=='valid':
                    valid[context['prefix_sha256']] = context['prefix_nt']
                else:
                    unresolved.append(context['status'])
            status = ('no_exact_source_match' if not row['canonical_matches'] else
                      'multiple_distinct_valid_prefixes' if len(valid)>1 else
                      'valid_context_with_unresolved_matches' if valid and unresolved else
                      'unique_valid_context' if valid else
                      'missing_leader_coordinates' if set(unresolved)=={'missing_leader_coordinates'} else
                      'invalid_context')
            coverage.append(dict(locus=locus,native_allele=name,native_sequence_sha256=row['sequence_sha256'],
                source_match_count=len(row['canonical_matches']),valid_distinct_prefix_count=len(valid),
                unresolved_match_count=len(unresolved),status=status))
        summaries[locus] = dict(native_v_alleles=len(legacy),public_v_alleles=len(alleles),
            native_context_status=dict(Counter(r['status'] for r in coverage if r['locus']==locus)),
            public_context_status=dict(Counter(c['status'] for c in contexts.values())))
    output.mkdir(parents=True)
    def table(name, rows):
        fields = list(dict.fromkeys(k for r in rows for k in r))
        with (output/name).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
    table('native_v_coverage.tsv',coverage)
    table('exact_source_matches.tsv',matches)
    table('public_v_contexts.tsv',raw_rows)
    summary=dict(purpose='raw model-parameter capability audit; not simulated cohort or runtime cache',
                 source_manifest=specification,installed_reference_sha256=identity['genairr_reference_sha256'],loci=summaries)
    (output/'audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summaries


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-dir',type=Path)
    p.add_argument('--output',type=Path)
    p.add_argument('--self-check',action='store_true')
    a=p.parse_args()
    if a.self_check:
        self_check();print('fixed transcript-coordinate checks passed')
    else:
        if a.reference_dir is None or a.output is None:p.error('--reference-dir and --output required')
        print(json.dumps(audit(a.reference_dir,a.output),indent=2))
