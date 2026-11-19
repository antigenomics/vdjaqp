"""Export observed preimputation contigs using the existing native adapter.

This stage consumes completed assembly observations only; it never generates a
repertoire or reads simulation truth. RNA mature boundaries require an explicit
public context catalogue, not generating-chain offsets.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def export_observations(root, geometry, output, *, cpus='all', boundary_contexts=None):
    from vdjaqp.resources import resolve_cpu_budget
    budget = resolve_cpu_budget(cpus)
    import polars as pl
    from vdjaqp.contigs import adapt_contigs

    root, output = Path(root), Path(output)
    if output.exists():
        raise FileExistsError('output already exists; choose a fresh directory')
    if geometry not in ('amplicon', 'rnaseq'):
        raise ValueError('unknown observation geometry')
    if boundary_contexts is not None and geometry != 'rnaseq':
        raise ValueError('boundary contexts apply only to rnaseq')
    sources = dict(raw=root/'assembly/raw_support.json',
                   annotations=root/'assembly/arda_contigs.json',
                   support=root/f'assembly/read_support/{geometry}_support.json')
    if boundary_contexts is not None:
        sources['boundary_contexts'] = Path(boundary_contexts)
    receipt_sources = {}
    loaded = {}
    for name, path in sources.items():
        data = path.read_bytes()
        receipt_sources[name] = dict(path=str(path.resolve()), bytes=len(data),
                                    sha256=hashlib.sha256(data).hexdigest())
        loaded[name] = json.loads(data)
    prefix = geometry + '_'
    raw = [r for r in loaded['raw'] if r['sequence_id'].startswith(prefix)]
    names = {r['sequence_id'] for r in raw}
    annotations = [r for r in loaded['annotations'] if r['sequence_id'] in names]
    supports = {}
    for row in loaded['support']:
        name = prefix + row['contig_id']
        if name in supports:
            raise ValueError('duplicate support contig_id')
        supports[name] = row
    if not raw:
        raise ValueError('no observed contigs for selected geometry')
    started = time.perf_counter()
    boundaries = None
    if boundary_contexts is not None:
        from vdjaqp.transcript_boundary import search_transcript_boundaries
        output.mkdir(parents=True, exist_ok=False)
        boundaries = search_transcript_boundaries(
            raw, annotations, supports, loaded['boundary_contexts'], threads=budget,
            alignment_output=output/'observed_boundary_alignments.tsv.gz')
        boundaries.write_parquet(output/'observed_transcript_boundaries.parquet')
    adapted = adapt_contigs(raw, annotations, supports,
                            **({'boundaries': boundaries} if boundaries is not None else {}))
    if boundary_contexts is None:
        output.mkdir(parents=True, exist_ok=False)
    counts = {}
    for name, frame in adapted.items():
        # Membership comes from observation provenance, never ID or assay truth.
        frame = frame.with_columns(pl.lit(True).alias('background_member'))
        frame.write_parquet(output/(name + '.parquet'))
        counts[name] = frame.height
    for name, rows in [('raw', raw), ('annotations', annotations), ('support', supports)]:
        (output/(name + '.json')).write_text(json.dumps(rows) + '\n')
    code = [Path(__file__), Path(sys.modules[adapt_contigs.__module__].__file__)]
    if boundaries is not None:
        code.append(Path(sys.modules[search_transcript_boundaries.__module__].__file__))
    receipt = dict(geometry=geometry, cpu_requested=str(cpus), cpu_resolved=budget,
                   sources=receipt_sources, counts=counts,
                   observed_locus_counts=dict(adapted['catalog'].group_by('locus').len().iter_rows()),
                   boundary_contexts_used=boundary_contexts is not None,
                   sequence_origin='native_preimputation_contig_retained_primary_support',
                   background_membership='observed_source_no_inserted_references',
                   code_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in code},
                   seconds=time.perf_counter()-started)
    (output/'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    return receipt


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--geometry', choices=('amplicon', 'rnaseq'), required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpus', default='all')
    parser.add_argument('--boundary-contexts', type=Path)
    args = parser.parse_args(argv)
    export_observations(args.root, args.geometry, args.output, cpus=args.cpus,
                        boundary_contexts=args.boundary_contexts)


if __name__ == '__main__':
    main()
