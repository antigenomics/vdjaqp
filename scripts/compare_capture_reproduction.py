"""Verify fresh capture outputs against all surviving accepted compact tables."""
import argparse
import hashlib
from pathlib import Path
import polars as pl

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('fresh',type=Path);p.add_argument('accepted',type=Path);p.add_argument('output',type=Path)
a=p.parse_args()
if a.output.exists():p.error('fresh output required')
cases=sorted(x for x in a.fresh.iterdir() if x.is_dir())
if {x.name for x in cases}!={'gel','solution',*[f'individual_{i:02d}' for i in range(1,9)]}:raise ValueError('complete context set required')
keys=dict(record_capture=['view','region','population'],gene_capture=['view','population','segment'],genes=['view','population','segment','gene_call'],junction_capture=['view','population'],mapping_ambiguity_distribution=['view','compatible_candidate_records','compatible_anchored_junction_strings'])
verdicts=[]
for name,group in keys.items():
    frames=[]
    for case in cases:
        source='junctions' if name=='junction_capture' else 'mapping_ambiguity' if name=='mapping_ambiguity_distribution' else name
        frame=pl.read_parquet(case/(source+'.parquet'))
        if name=='junction_capture':
            frame=frame.group_by(group).agg(pl.len().alias('unique_anchored_junction_strings'),pl.col('any_partial_capture').sum().alias('any_junction_capture'),(pl.col('any_partial_capture')&~pl.col('complete_capture')).sum().alias('partial_only_junction_capture'),pl.col('complete_capture').sum().alias('complete_junction_capture'),pl.col('member_records').sum()).with_columns((pl.col('any_junction_capture')/pl.col('unique_anchored_junction_strings')).alias('any_capture_fraction'),(pl.col('complete_junction_capture')/pl.col('unique_anchored_junction_strings')).alias('complete_capture_fraction'))
        elif name=='mapping_ambiguity_distribution':
            frame=frame.group_by(group).len(name='peptide_observation_rows')
        frames.append(frame.with_columns(pl.lit(case.name).alias('case')))
    actual=pl.concat(frames)
    path=a.accepted/(name+'.tsv.gz')
    baseline=pl.read_csv(path,separator='\t',null_values=['NA','null'])
    for col in baseline.columns:
        if actual.schema[col]==pl.Null:
            if baseline[col].null_count()!=baseline.height:raise ValueError('lost non-null measurement')
            baseline=baseline.with_columns(pl.lit(None).alias(col))
        else:baseline=baseline.with_columns(pl.col(col).cast(actual.schema[col],strict=True))
    sort=['case',*group]
    actual=actual.select(baseline.columns).sort(sort);baseline=baseline.sort(sort)
    if not actual.equals(baseline):raise ValueError(f'changed accepted capture values: {name}')
    verdicts.append(dict(table=name,rows=actual.height,columns=actual.width,
        exact_values_equal=True,missing_values_preserved=True,
        accepted_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        canonical_tsv_sha256=hashlib.sha256(actual.write_csv(separator='\t').encode()).hexdigest()))
a.output.mkdir()
pl.DataFrame(verdicts).write_csv(a.output/'capture_comparison.tsv',separator='\t')
print(f'{len(verdicts)} compact tables exactly reproduced across ten contexts')
