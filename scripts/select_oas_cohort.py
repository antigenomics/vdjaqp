"""Select the authorized OAS protein cohort afresh; retain input provenance and attrition."""
import argparse
import hashlib
import json
from pathlib import Path
import polars as pl

p = argparse.ArgumentParser()
p.add_argument('--input', type=Path, required=True)
p.add_argument('--influenza-metadata', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
if a.output.exists():
    raise SystemExit('Output must be new; no resume or computed-input reuse.')
meta = pl.read_csv(a.influenza_metadata, separator='\t', infer_schema=False)
keys = (meta.filter((pl.col('study') == 'Tipton_2015') &
                   (pl.col('vaccine') == 'Flu') &
                   pl.col('subject').str.starts_with('Subject-FLU'))
        .select((pl.col('study') + '|' + pl.col('subject')).alias('donor'))
        .unique()['donor'].to_list())
assert len(keys) == 4
f = pl.scan_parquet(a.input).with_columns(
    pl.when(pl.col('disease') == 'SARS-COV-2').then(pl.lit('covid_disease'))
    .when(pl.col('donor').is_in(keys)).then(pl.lit('influenza_vaccination_associated'))
    .otherwise(pl.lit(None)).alias('cohort'))
f = f.filter(pl.col('cohort').is_not_null()).collect()
checks = {
    'canonical_amino_acids': pl.col('sequence').str.contains('^[ACDEFGHIKLMNPQRSTVWY]+$'),
    'heavy_vj_calls': pl.col('v_gene').str.starts_with('IGHV') & pl.col('j_gene').str.starts_with('IGHJ'),
    'consistent_lengths': (pl.col('sequence').str.len_chars() == pl.col('length')) & (pl.col('region_ids').list.len() == pl.col('length')),
    'region_labels_present': pl.col('region_ids').list.eval(pl.element().is_not_null() & (pl.element() >= 0)).list.all(),
}
counts = []
for name, rule in [('metadata_eligible', pl.lit(True)), *checks.items()]:
    f = f.filter(rule.fill_null(False))
    counts.extend(dict(stage=name, **row) for row in f.group_by('cohort').len().rename({'len':'records'}).sort('cohort').to_dicts())
assert f.height and f['clone_id'].n_unique() == f.height
f = f.sort('cohort', 'study', 'donor', 'clone_id')
a.output.mkdir(parents=True)
f.write_parquet(a.output / 'records.parquet', compression='zstd')
pl.DataFrame(counts).write_csv(a.output / 'attrition.tsv', separator='\t')
f.group_by('cohort','study','donor').len().rename({'len':'records'}).sort('cohort','study','donor').write_csv(a.output/'populations.tsv',separator='\t')
def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()
receipt = dict(input_sha256=sha(a.input), metadata_sha256=sha(a.influenza_metadata),
               script_sha256=sha(Path(__file__)), polars_version=pl.__version__,
               output_sha256=sha(a.output/'records.parquet'), records=f.height,
               unique_amino_acid_sequences=f['sequence'].n_unique(),
               checks=list(checks), influenza_donor_keys=sorted(keys),
               interpretation='New sequence-integrity-qualified protein cohort; all qualifying source records retained, no random sampling or deduplication.',
               limitations=['Not the historical 473000 subset.', 'No nucleotide sequences, productive flag or raw-read quality available.', 'Region-label semantics and complete variable-domain boundaries are not certified.', 'Disease/vaccination association does not establish antigen specificity.', 'Influenza selection is restricted to four explicitly labelled Tipton donors; timepoints cannot be resolved in this export.'])
(a.output/'provenance.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'records':f.height,'counts':counts},indent=2))
