"""Rebuild the authorized OAS selection and annotate observed proteins with native ARDA.

Run as ``python -m scripts.annotate_oas_cohort`` in the pinned project environment.
This fresh-run driver preserves unresolved annotations and never reconstructs DNA.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

from vdjaqp.resources import resolve_cpu_budget


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--influenza-metadata', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpus', default='all')
    args = parser.parse_args()
    cpus = resolve_cpu_budget(args.cpus)
    import arda
    import polars as pl
    import pyarrow as pa
    import pyarrow.csv as csv
    import pyarrow.parquet as pq

    if args.output.exists():
        parser.error('Output must be new; no resume or computed-input reuse.')
    if arda.__version__ != '2.38.7':
        raise ValueError('This declared analysis requires pinned ARDA 2.38.7.')
    started = time.monotonic()
    source = Path(__file__).resolve()
    selector = source.with_name('select_oas_cohort.py')
    pins = {str(p.resolve()): sha(p) for p in (args.input, args.influenza_metadata, source, selector,
                                             source.parents[1] / 'uv.lock', source.parents[1] / 'pyproject.toml')}
    args.output.mkdir(parents=True)
    selection = args.output / 'selection'
    subprocess.run([sys.executable, str(selector), '--input', str(args.input),
                    '--influenza-metadata', str(args.influenza_metadata),
                    '--output', str(selection)], check=True)
    proteins = pl.read_parquet(selection / 'records.parquet', columns=['clone_id', 'sequence'])
    fasta = args.output / 'proteins.fasta'
    proteins.select(pl.concat_str(pl.lit('>'), 'clone_id', pl.lit('\n'), 'sequence')).write_csv(
        fasta, include_header=False, quote_style='never')
    reference = arda.reference_index_identity('human', seqtype='aa')
    tsv = args.output / 'annotations.airr.tsv'
    command = [str(Path(sys.executable).parent / 'arda'), '--log-file', str(args.output / 'arda.log'),
               'annotate', '--input', str(fasta), '--output', str(tsv), '--seqtype', 'aa',
               '--organism', 'human', '--threads', str(cpus)]
    subprocess.run(command, check=True)
    if reference != arda.reference_index_identity('human', seqtype='aa'):
        raise ValueError('Native reference changed during annotation.')
    # Preserve native AIRR text semantics and empty fields, using bounded Arrow batches.
    with tsv.open() as stream:
        names = stream.readline().rstrip('\n').split('\t')
    annotated = args.output / 'annotations.parquet'
    reader = csv.open_csv(tsv, read_options=csv.ReadOptions(block_size=16 * 1024 * 1024),
                          parse_options=csv.ParseOptions(delimiter='\t'),
                          convert_options=csv.ConvertOptions(column_types={n: pa.string() for n in names}))
    with pq.ParquetWriter(annotated, reader.schema, compression='zstd') as writer:
        for batch in reader:
            writer.write_batch(batch)
    rows = pl.read_parquet(annotated, columns=['sequence_id', 'sequence'])
    if rows.height != proteins.height or rows['sequence_id'].n_unique() != proteins.height:
        raise ValueError('Native annotation lost or duplicated input records.')
    compared = proteins.join(rows, left_on='clone_id', right_on='sequence_id', how='left', validate='1:1')
    if compared.filter(pl.col('sequence_right').is_null() | (pl.col('sequence') != pl.col('sequence_right'))).height:
        raise ValueError('Native annotation changed or omitted observed protein sequences.')
    if any(sha(Path(p)) != expected for p, expected in pins.items()):
        raise ValueError('Input or source changed during the fresh run.')
    receipt = dict(status='complete', records=proteins.height, cpus=cpus,
                   seconds=time.monotonic() - started, arda_version=arda.__version__,
                   source_input_sha256=pins, reference=reference, command=command,
                   annotations_sha256=sha(annotated), fasta_sha256=sha(fasta), airr_sha256=sha(tsv),
                   semantics='Observed amino-acid input; native inferred regions and D calls. No reconstructed nucleotide sequence.',
                   retention='Original input and final selection/annotations retained; redundant FASTA and AIRR text removed after exact record/sequence verification.')
    fasta.unlink()
    tsv.unlink()
    (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': 'complete', 'records': proteins.height, 'seconds': receipt['seconds']}))


if __name__ == '__main__':
    main()
