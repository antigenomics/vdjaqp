"""Build requested donor-context model parameters from verified public references.

Use --prepare in the pinned GenAIRR environment, then --align in the ARDA
environment. This is a reference capability audit, never a repertoire draw/cache.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))


def prepare(raw_dir, output):
    import pickle
    from audit_transcript_context import load_verified_sources, extract_context, normalized
    from genairr_gene_map import exact_v_map
    manifest, identity, raw, pickles = load_verified_sources(raw_dir)
    native, sources = [], []
    for locus, alleles in raw.items():
        config = pickle.loads(pickles[locus])
        legacy = [(a.name, a.gene, normalized(a.ungapped_seq))
                  for group in config.v_alleles.values() for a in group]
        canonical = [(a['label'], a['label'].split('*')[0], normalized(a['coding_sequence']))
                     for a in alleles]
        mapping = exact_v_map(legacy, canonical)
        for name, gene, sequence in legacy:
            native.append(dict(locus=locus, native_allele=name, gene=gene, core_nt=sequence,
                               exact_sources=[f'{locus}:{sid}' for sid in mapping[name]['canonical_alleles']]))
        for allele in alleles:
            sources.append(dict(source_id=f"{locus}:{allele['label']}", locus=locus,
                                source_allele=allele['label'],
                                core_nt=normalized(allele['coding_sequence']), **extract_context(allele)))
    output = Path(output)
    if output.exists():
        raise ValueError('model source output must be fresh')
    output.write_text(json.dumps(dict(purpose='verified raw-reference model parameters; not cohort intermediates',
        source_manifest=manifest, native_reference_sha256=identity['genairr_reference_sha256'],
        native=native, sources=sources), indent=2) + '\n')


def align(source_path, output, threads):
    from arda import mmseqs
    from vdjaqp.transcript_context import select_contexts
    from vdjaqp.resources import resolve_cpu_budget
    threads = resolve_cpu_budget(threads)
    output = Path(output)
    if output.exists():
        raise ValueError('capability audit output must be fresh')
    source_bytes = Path(source_path).read_bytes()
    model = json.loads(source_bytes)
    native, sources = model['native'], model['sources']
    first, provenance = select_contexts(native, sources)
    strict = {r['native_allele'] for r in first if r['method'] == 'STRICT_EXACT'}
    donors = {p['source_id'] for p in provenance if p['native_allele'] in strict}
    missing = {r['native_allele'] for r in first if r['method'] == 'UNSUPPORTED'}
    output.mkdir(parents=True)
    query, target = output/'native_v.fasta', output/'donor_v.fasta'
    query.write_text(''.join(f">{r['native_allele']}\n{r['core_nt']}\n" for r in native if r['native_allele'] in missing))
    target.write_text(''.join(f">{r['source_id']}\n{r['core_nt']}\n" for r in sources if r['source_id'] in donors))
    if len(donors) > 300:
        raise ValueError('donor count exceeds unchanged native search candidate cap')
    if __import__('os').environ.get('ARDA_MMSEQS_SEARCH_OPTS', '').strip():
        raise ValueError('unsupported MMseqs override')
    qdb = mmseqs.createdb(query, output/'queries', dbtype=2)
    tdb = mmseqs.createdb(target, output/'targets', dbtype=2)
    result = mmseqs.search(qdb, tdb, output/'results', output/'work', search_type=3,
        sensitivity=7.0, evalue=1e-3, max_seqs=300, threads=threads, kmer=13, extra=['--strand', '1'])
    table = mmseqs.convertalis(qdb, tdb, result, output/'alignments.tsv', search_type=3, threads=threads)
    ns = {r['native_allele']: r for r in native}
    ss = {r['source_id']: r for r in sources}
    hits, rejected_cross_locus = [], 0
    with table.open() as stream:
        for h in csv.DictReader(stream, fieldnames=mmseqs.DEFAULT_FORMAT_OUTPUT.split(','), delimiter='\t'):
            if int(h['qstart']) > int(h['qend']) or int(h['tstart']) > int(h['tend']):
                raise ValueError('unexpected reverse alignment')
            if ns[h['query']]['locus'] != ss[h['target']]['locus']:
                rejected_cross_locus += 1
                continue
            hits.append(dict(native_allele=h['query'], source_id=h['target'], bits=float(h['bits'])))
    selections, candidates = select_contexts(native, sources, hits)
    def write(name, rows):
        with (output/name).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader(); writer.writerows(rows)
    write('contexts.tsv', selections)
    write('source_candidates.tsv', candidates)
    selected = {r['native_allele']: r for r in selections}
    boundary_contexts = [dict(target=r['native_allele'],
        sequence=selected[r['native_allele']]['prefix_nt'] + r['core_nt'],
        mature_start=len(selected[r['native_allele']]['prefix_nt']) + 1)
        for r in native if selected[r['native_allele']]['prefix_nt'] is not None]
    (output/'boundary_contexts.json').write_text(json.dumps(boundary_contexts, indent=2)+'\n')
    summary = dict(purpose='source-only leader imputation capability; no cohort or draw',
        model_source_sha256=hashlib.sha256(source_bytes).hexdigest(), native_alleles=len(native),
        initial_methods=dict(Counter(r['method'] for r in first)),
        final_methods=dict(Counter(r['method'] for r in selections)),
        confidence=dict(Counter(r['confidence'] for r in selections)),
        donor_source_count=len(donors), same_locus_alignment_count=len(hits),
        rejected_cross_locus_alignment_count=rejected_cross_locus, threads=threads,
        search=dict(sensitivity=7.0, evalue=1e-3, max_seqs=300, kmer=13, strand='forward'))
    (output/'audit.json').write_text(json.dumps(summary, indent=2)+'\n')
    # Search DBs and native work are task scratch, not retained computational caches.
    import shutil
    for path in output.iterdir():
        if path.name not in {'native_v.fasta', 'donor_v.fasta', 'alignments.tsv', 'contexts.tsv', 'source_candidates.tsv', 'audit.json', 'boundary_contexts.json'}:
            if path.is_dir(): shutil.rmtree(path)
            else: path.unlink()
    print(json.dumps(summary, indent=2))
    if summary['final_methods'].get('UNSUPPORTED'):
        raise RuntimeError('native alleles without donor hits; inspect retained source metadata')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    modes = p.add_mutually_exclusive_group(required=True)
    modes.add_argument('--prepare', type=Path, help='verified raw public reference directory')
    modes.add_argument('--align', type=Path, help='prepared verified model source metadata')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--threads', type=int, default=16)
    a = p.parse_args()
    if a.prepare is not None:
        prepare(a.prepare, a.output)
    else:
        if a.threads < 1: p.error('threads must be positive')
        align(a.align, a.output, a.threads)
