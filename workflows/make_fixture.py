"""Create a tiny technical acceptance fixture, never a scientific benchmark."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import polars as pl
from vdjaqp.errormodel import mass_prior

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("output", type=Path)
parser.add_argument("--mock-ms", action="store_true", help="Use the bundled 10,000-call MS fixture")
args = parser.parse_args()
out = args.output.resolve()
out.mkdir(parents=True, exist_ok=False)
pl.DataFrame({'clone_id':['c0'], 'seq':['ACDEFGHKLMNP'],
    'cdr3_seq':['ACDEFGHKLMNP'], 'region_ids':[[6]*12],
    'cdr3_start':[0], 'cdr3_end':[12], 'background_member':[True],
    'v_name':['IGHV1-2'], 'j_name':['IGHJ4'], 'c_name':['IGHG1']}).write_parquet(out/'candidates.parquet')
calls = (
    pl.read_csv(Path(__file__).resolve().parents[1] / 'examples/inputs/mock-ms.tsv.gz',
                separator='\t')
    if args.mock_ms else pl.DataFrame({'peptide_id': ['run/100/2/50/10/1/ACDEFG'],
                                     'seq': ['ACDEFG']})
)
calls.with_columns(
    pl.lit([float('nan')] * 6).alias('local_confidence'),
    pl.lit('db').alias('search_type'),
    pl.lit(float('nan')).alias('alc'),
    pl.lit(1.).alias('alc_rank'),
    pl.lit(1, dtype=pl.Int64).alias('n_reads'),
).write_parquet(out / 'peptides.parquet')
(out/'model.json').write_text(json.dumps(dict(conf_x=[0.,1.],conf_e=[.1,.01],
    substitution=mass_prior().tolist(), e_db=.1)))
with (out / 'samples.csv').open('w', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['sample', 'candidates', 'peptides', 'producer', 'model'])
    writer.writerow(['technical', out / 'candidates.parquet', out / 'peptides.parquet',
                     'technical fixture', out / 'model.json'])
