"""Compare complete fresh real metrics to surviving original/corrected source tables."""
import argparse
import hashlib
import json
from pathlib import Path
import polars as pl

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('fresh',type=Path);p.add_argument('accepted',type=Path);p.add_argument('output',type=Path)
a=p.parse_args()
if a.output.exists():p.error('fresh output required')
f=pl.read_csv(a.fresh,separator='\t');b=pl.read_csv(a.accepted,separator='\t')
cases={'gel','solution',*[f'individual_{i:02d}' for i in range(1,9)]}
methods={'alignment_union','exclusive_spectrum'}
columns='n_candidates n_positives n_benchmark_negatives average_precision roc_auc threshold_confirmed_nats true_positive false_positive false_negative true_negative precision_at20 recall_at20 f1_at20 n_clonotypes median_rank worst_rank top_50 top_500'.split()
counts={'true_positive','false_positive','false_negative','true_negative','top_50','top_500'}
for frame in (f,b):
    if not set(columns)<=set(frame.columns):raise ValueError('all 18 canonical metric columns required')
    if frame.select(pl.col(columns).null_count()).sum_horizontal().item() or not all(frame.select(pl.col(c).cast(pl.Float64).is_finite().all() for c in columns).row(0)):
        raise ValueError('all canonical metrics must be finite/non-null')
    for c in columns:
        if c.startswith('n_') or c in counts or c=='worst_rank':
            if frame.filter(pl.col(c)!=pl.col(c).floor()).height:raise ValueError('counts and ordinal ranks must be exact integers')
if (f.height!=20 or set(f['case'])!=cases or set(f['method'])!=methods
        or f.select('case','method').unique().height!=20):
    raise ValueError('exact ten-context/two-method keys required')
expected_methods={'original arithmetic alignment_union','corrected arithmetic alignment_union','corrected arithmetic exclusive_spectrum'}
if b.height!=30 or set(b['case'])!=cases or set(b['method'])!=expected_methods or b.select('case','method').unique().height!=30:
    raise ValueError('complete original/corrected accepted endpoint keys required')
b=b.with_columns(pl.when(pl.col('method').str.ends_with('exclusive_spectrum')).then(pl.lit('exclusive_spectrum')).otherwise(pl.lit('alignment_union')).alias('current_method'))
joined=b.join(f.rename({'method':'current_method'}),on=['case','current_method'],suffix='_fresh',validate='m:1')
if joined.height!=b.height:raise ValueError('accepted populations not matched')
rows=[]
for row in joined.iter_rows(named=True):
    for metric in columns:
        old,new=row[metric],row[metric+'_fresh']
        count=metric.startswith('n_') or metric in counts
        rows.append(dict(case=row['case'],accepted_method=row['method'],current_method=row['current_method'],
            quantity=metric,accepted_value=str(int(old)) if count else str(old),
            current_value=str(int(new)) if count else str(new),delta=new-old,
            units='rank' if 'rank' in metric else 'count' if count else 'nats' if metric=='threshold_confirmed_nats' else 'fraction',
            n_candidates=row['n_candidates'],n_positives=row['n_positives']))
a.output.mkdir(parents=True)
out=pl.DataFrame(rows);out.write_csv(a.output/'matched_metric_deltas.tsv',separator='\t')
corrected=out.filter(pl.col('accepted_method').str.starts_with('corrected'))
report=dict(contexts=10,current_method_contexts=f.height,matched_metric_rows=out.height,
    corrected_metric_rows=corrected.height,corrected_max_abs_delta=corrected['delta'].abs().max(),
    corrected_counts_exact=corrected.filter((pl.col('units')=='count') & (pl.col('delta')!=0)).height==0,
    corrected_all_equal_at1e12=corrected.filter(pl.col('delta').abs()>1e-12).height==0,
    original_baseline_preserved=True,inputs_sha256={n:hashlib.sha256(path.read_bytes()).hexdigest() for n,path in [('fresh',a.fresh),('accepted',a.accepted)]})
(a.output/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
