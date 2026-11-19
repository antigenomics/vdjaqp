"""Report native parent-clone assignment capture from final VDJ-AQP deliverables.

Collapse compatible variants before joining measured native parent counts, so
one clone's counts are never multiplied over its candidate variants. Unassessed
clones remain unknown. This is a descriptive source-assignment readout, not
protein abundance or accepted identity, and introduces no new threshold.
"""
import argparse,json,hashlib
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--project-root',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--cpus',type=int)
a=p.parse_args()
import vdjaqp
budget=vdjaqp.resolve_cpu_budget(a.cpus)
import polars as pl
if a.output.exists():raise FileExistsError(a.output)
summary=[];receipts=[]
cases=['gel','solution']+[f'individual_{i:02d}' for i in range(1,9)]
for cohort in ['S1','pooled']:
 inputs=a.project_root/'protected_inputs'/cohort/'inputs'
 ledger=pl.read_parquet(inputs/'native_clone_ledger.parquet').filter(pl.col('locus')=='IGH')
 mapping=pl.read_parquet(inputs/'heavy_variants.parquet').select('clone_id','clone_row')
 assert mapping['clone_id'].n_unique()==mapping.height and set(mapping['clone_row'])<=set(ledger['clone_row'])
 counts=ledger.select('clone_row',*[pl.col(k).cast(pl.UInt64,strict=True) for k in ['duplicate_count','umi_count']])
 assert all(counts[k].null_count()==0 for k in ['duplicate_count','umi_count'])
 for case in cases:
  source=a.project_root/'protected_context_analysis'/cohort/case/'analysis'
  provenance=json.loads((source/'provenance.json').read_text());assert provenance['status']=='complete'
  path=source/case/'records.parquet'
  assert hashlib.sha256(path.read_bytes()).hexdigest()==provenance['artifact_sha256'][f'{case}/records.parquet']
  d=pl.scan_parquet(path).filter((pl.col('region')=='protein')&~pl.col('is_reference')).select('view','clone_id','covered_aa','defined_aa').collect()
  assert d.filter((pl.col('defined_aa')<=0)|(pl.col('covered_aa')<0)|(pl.col('covered_aa')>pl.col('defined_aa'))).height==0
  views=d['view'].unique().to_list();assert len(views)==3
  for view in views:
   part=d.filter(pl.col('view')==view)
   assert part['clone_id'].n_unique()==mapping.height and part.height==mapping.height and set(part['clone_id'])==set(mapping['clone_id'])
  collapsed=d.join(mapping,on='clone_id',validate='m:1').group_by('view','clone_row').agg(
    pl.len().alias('assessed_variant_records'),(pl.col('covered_aa')>0).any().alias('any_variant_covered'))
  grid=counts.join(pl.DataFrame({'view':views}),how='cross').join(collapsed,on=['view','clone_row'],how='left',validate='1:1')
  grid=grid.with_columns(pl.col('assessed_variant_records').is_not_null().alias('assessed'))
  for view in views:
   block=grid.filter(pl.col('view')==view)
   for column,unit in [('duplicate_count','consensus_assignments'),('umi_count','distinct_UMI_assignments')]:
    total=sum(map(int,block[column]));assessed=sum(map(int,block.filter(pl.col('assessed'))[column]));captured=sum(map(int,block.filter(pl.col('any_variant_covered').fill_null(False))[column]))
    assert 0<=captured<=assessed<=total<2**63
    summary.append(dict(cohort=cohort,case=case,view=view,count_unit=unit,native_heavy_clones=block.height,
      assessed_clones=int(block['assessed'].sum()),unassessed_clones=int((~block['assessed']).sum()),
      captured_clones=int(block['any_variant_covered'].fill_null(False).sum()),all_source_assignments=total,
      assessed_source_assignments=assessed,unassessed_source_assignments=total-assessed,captured_source_assignments=captured,
      captured_fraction_of_assessed=captured/assessed if assessed else None,observed_capture_lower_bound_of_all=captured/total if total else None))
  dest=a.output/cohort;dest.mkdir(parents=True,exist_ok=True)
  grid.write_parquet(dest/f'{case}.parquet',compression='zstd')
  receipts.append(dict(cohort=cohort,case=case,record_capture_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
pl.DataFrame(summary).write_csv(a.output/'native_assignment_capture.tsv',separator='\t')
(a.output/'manifest.json').write_text(json.dumps(dict(status='complete',contexts=20,cpus=budget,
 unit='native corrected heavy-clone source assignments; not global molecules, transcript abundance or protein amounts',
 unknown='parents without tested observed protein variants retain null capture and are reported separately',
 evidence='existing final capture views; compatibility/positive LLR/conditional 20-nat view do not establish unique identity',
 receipts=receipts,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
print('Complete: 20 contexts, native parent counts collapsed once per view; unassessed parents explicit')
