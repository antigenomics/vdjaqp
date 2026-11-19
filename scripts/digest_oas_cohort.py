"""Aggregate a freshly completed OAS protein annotation using qualified native OpenMS.

Final outputs are integer peptide-length counts and internal-bond cleavage maps.
Unknown inferred regions remain in native-position and peptide-length summaries.
"""
import argparse
import json
from pathlib import Path
import subprocess
import time

from vdjaqp.resources import resolve_cpu_budget
from scripts.annotate_oas_cohort import sha

REGIONS = ('fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'fwr4')
ENZYMES = ('Trypsin', 'Chymotrypsin', 'Asp-N')


def prepare(proteins, annotations):
    """Validate native coordinates against observed substrings before changing conventions."""
    import polars as pl
    if annotations.height != proteins.height or annotations['sequence_id'].n_unique() != proteins.height:
        raise ValueError('Annotation records are missing or duplicated.')
    frame = proteins.join(annotations, left_on='clone_id', right_on='sequence_id', how='left', validate='1:1')
    if frame.filter(pl.col('sequence_right').is_null() | (pl.col('sequence') != pl.col('sequence_right'))).height:
        raise ValueError('Observed protein/annotation identity mismatch.')
    frame = frame.drop('sequence_right').with_columns(pl.col('sequence').str.len_chars().alias('length'))
    previous = pl.lit(0)
    for region in REGIONS:
        start, end, aa = (region+s for s in ('_start', '_end', '_aa'))
        frame = frame.with_columns(pl.col(start, end).cast(pl.String).replace('', None).cast(pl.Int64), pl.col(aa).fill_null(''))
        present = pl.col(start).is_not_null() & pl.col(end).is_not_null()
        valid = ((pl.col(start) >= 1) & (pl.col(end) >= pl.col(start)) & (pl.col(end) <= pl.col('length'))
                 & (pl.col(start) > previous)
                 & (pl.col('sequence').str.slice(pl.col(start)-1, pl.col(end)-pl.col(start)+1) == pl.col(aa)))
        if frame.filter((present & ~valid.fill_null(False)) |
                        (pl.col(start).is_null() != pl.col(end).is_null()) |
                        (~present & (pl.col(aa) != ''))).height:
            raise ValueError(f'Inconsistent native {region} coordinates or substring.')
        previous = pl.coalesce(pl.col(end), previous)
    frame = frame.with_columns(pl.sum_horizontal([pl.col(r+'_end').is_not_null().cast(pl.Int64) for r in REGIONS]).alias('defined_regions'))
    quality = frame.group_by('cohort', 'defined_regions').len(name='protein_records').sort('cohort', 'defined_regions')
    cohorts = frame.select('cohort').unique().sort('cohort').with_row_index('cohort_id')
    genes = frame.select(pl.col('v_call').fill_null('').alias('v_call')).unique().sort('v_call').with_row_index('gene_id')
    frame = frame.with_columns(pl.col('v_call').fill_null('')).join(cohorts, on='cohort', validate='m:1').join(genes, on='v_call', validate='m:1')
    fields = ['cohort_id', 'gene_id', 'sequence']
    for region in REGIONS:
        fields += [(pl.col(region+'_start')-1).fill_null(0).alias(region+'_start0'),
                   pl.col(region+'_end').fill_null(0).alias(region+'_end0')]
    return frame.select(fields), cohorts, genes, quality


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotations', type=Path, required=True, help='Completed fresh annotation directory')
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--native-receipt', type=Path, required=True)
    parser.add_argument('--openms-prefix', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpus', default='all')
    args = parser.parse_args()
    cpus = resolve_cpu_budget(args.cpus)
    import polars as pl
    if args.output.exists():
        parser.error('Output must be new; no resume.')
    started = time.monotonic()
    native = json.loads(args.native_receipt.read_text())
    source = Path(__file__).with_name('oas_digest.cpp')
    library = args.openms_prefix/'lib/libOpenMS.so'
    enzymes = args.openms_prefix/'share/OpenMS/CHEMISTRY/Enzymes.xml'
    for path, key in [(source, 'source_sha256'), (args.binary, 'binary_sha256'),
                      (library, 'library_sha256'), (enzymes, 'enzyme_definitions_sha256')]:
        if sha(path) != native[key]:
            raise ValueError(f'Qualified native artifact changed: {key}')
    receipt_path = args.annotations/'receipt.json'
    annotation_receipt = json.loads(receipt_path.read_text())
    annotations_path = args.annotations/'annotations.parquet'
    if annotation_receipt['status'] != 'complete' or sha(annotations_path) != annotation_receipt['annotations_sha256']:
        raise ValueError('Incomplete or changed annotation deliverable.')
    selection = args.annotations/'selection'
    selection_receipt = json.loads((selection/'provenance.json').read_text())
    proteins_path = selection/'records.parquet'
    # The selector pins the actual final table, independently of the annotation receipt.
    if sha(proteins_path) != selection_receipt['output_sha256']:
        raise ValueError('Selected protein table changed.')
    inputs = [receipt_path, annotations_path, proteins_path, args.native_receipt, args.binary,
              source, Path(__file__), library, enzymes, selection/'provenance.json']
    pins = {str(p.resolve()): sha(p) for p in inputs}
    proteins = pl.read_parquet(proteins_path, columns=['clone_id', 'cohort', 'sequence'])
    columns = ['sequence_id', 'sequence', 'v_call'] + [r+s for r in REGIONS for s in ('_start', '_end', '_aa')]
    annotations = pl.read_parquet(annotations_path, columns=columns)
    table, cohorts, genes, quality = prepare(proteins, annotations)
    del annotations, proteins
    totals = table.group_by('cohort_id').agg(pl.len().alias('protein_records'), pl.col('sequence').str.len_chars().sum().alias('residues'))
    args.output.mkdir(parents=True)
    cohorts.write_csv(args.output/'cohorts.tsv', separator='\t')
    genes.write_csv(args.output/'genes.tsv', separator='\t')
    quality.write_csv(args.output/'annotation_quality.tsv', separator='\t')
    temporary = args.output/'proteins.tsv'
    table.write_csv(temporary, separator='\t', include_header=False, quote_style='never')
    del table
    command = [str(args.binary.resolve()), str(temporary), str(args.output/'counts-'), str(cpus)]
    subprocess.run(command, check=True)
    outputs = {}
    for name in ('native', 'aligned', 'lengths', 'records'):
        tsv = args.output/f'counts-{name}.tsv'
        frame = pl.read_csv(tsv, separator='\t')
        if name in ('native', 'aligned'):
            if frame.filter(pl.any_horizontal([(pl.col(e) < 0) | (pl.col(e) > pl.col('available_bonds')) for e in ENZYMES])).height:
                raise ValueError('Cleavage counts exceed available internal bonds.')
            if name == 'native':
                bonds = frame.group_by('cohort_id').agg(pl.col('available_bonds').sum())
                if bonds.height != totals.height or bonds.join(totals, on='cohort_id').filter(
                        pl.col('available_bonds') != pl.col('residues')-pl.col('protein_records')).height:
                    raise ValueError('Native map omitted internal bonds.')
        elif name == 'records':
            if not frame.sort('cohort_id').equals(totals.select('cohort_id', 'protein_records').cast(frame.schema).sort('cohort_id')):
                raise ValueError('Native digestion lost protein records.')
        elif name == 'lengths':
            zero = frame.filter(pl.col('missed_cleavages') == 0).group_by('cohort_id', 'enzyme').agg(
                (pl.col('length_aa')*pl.col('peptide_occurrences')).sum().alias('digested_residues'))
            if zero.height != 3*totals.height or zero.join(totals, on='cohort_id').filter(pl.col('digested_residues') != pl.col('residues')).height:
                raise ValueError('Zero-missed fragments do not recover all observed residues.')
            eligible = frame.group_by('cohort_id', 'enzyme', 'missed_cleavages').agg(
                pl.col('peptide_occurrences').sum().alias('all_peptide_occurrences'),
                pl.col('peptide_occurrences').filter(pl.col('length_aa').is_between(5, 50)).sum().alias('length_eligible_occurrences'))
            eligible.sort('cohort_id', 'enzyme', 'missed_cleavages').write_csv(args.output/'length_eligibility.tsv', separator='\t')
        target = args.output/f'{name}.parquet'
        frame.write_parquet(target, compression='zstd')
        outputs[target.name] = sha(target)
        tsv.unlink()
    if any(sha(Path(p)) != expected for p, expected in pins.items()):
        raise ValueError('Source changed during digestion.')
    temporary.unlink()
    outputs.update({p.name: sha(p) for p in args.output.glob('*.tsv')})
    report = dict(status='complete', records=int(totals['protein_records'].sum()), cpus=cpus,
                  seconds=time.monotonic()-started, command=command, source_input_sha256=pins,
                  output_sha256=outputs, native=native,
                  semantics='Distinct positional occurrences; exact missed-cleavage strata 0,1,2,3; fully specific native digestion. Internal-bond denominators exclude terminal ends. Region coordinates are native inferences; each bond belongs to its preceding residue. Unresolved regions remain in native maps and length counts.',
                  display='Native protein length bins: [10*b,10*b+9]; internal bond bins: [5*b+1,5*b+5]. Aligned maps use left-aligned region-local bond offsets and lexically sorted native V calls; empty call retained.',
                  interpretation='Length eligibility (5–50 residues) is not MS detectability; no mass window or reconstructed nucleotide sequence.')
    (args.output/'receipt.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'complete', 'records': report['records'], 'seconds': report['seconds']}))


if __name__ == '__main__':
    main()
