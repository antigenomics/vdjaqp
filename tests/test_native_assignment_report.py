"""Parent counts are exact, counted once, and unavailable capture stays null."""
from pathlib import Path
import hashlib,json,subprocess,sys
import polars as pl
def test_native_assignment_counts_are_not_multiplied_by_variants(tmp_path):
 root=tmp_path
 for cohort in ['S1','pooled']:
  inp=root/'protected_inputs'/cohort/'inputs';inp.mkdir(parents=True)
  pl.DataFrame(dict(clone_row=[0,1],locus=['IGH','IGH'],duplicate_count=[str(2**53+1),'7'],umi_count=['3','5'])).write_parquet(inp/'native_clone_ledger.parquet')
  pl.DataFrame(dict(clone_id=['a','b'],clone_row=[0,0])).write_parquet(inp/'heavy_variants.parquet')
  for case in ['gel','solution']+[f'individual_{i:02d}' for i in range(1,9)]:
   path=root/'protected_context_analysis'/cohort/case/'analysis';(path/case).mkdir(parents=True)
   rows=[dict(view=v,clone_id=c,region='protein',is_reference=False,covered_aa=1,defined_aa=3) for v in ['all_compatible_alignments','positive_llr_alignments','identified_candidates_at20'] for c in ['a','b']]
   data=path/case/'records.parquet';pl.DataFrame(rows).write_parquet(data)
   (path/'provenance.json').write_text(json.dumps(dict(status='complete',artifact_sha256={f'{case}/records.parquet':hashlib.sha256(data.read_bytes()).hexdigest()})))
 out=root/'report'
 subprocess.run([sys.executable,str(Path(__file__).resolve().parents[1]/'scripts/report_native_assignments.py'),'--project-root',str(root),'--output',str(out)],check=True)
 data=pl.read_csv(out/'native_assignment_capture.tsv',separator='\t').filter(pl.col('count_unit')=='consensus_assignments')
 assert (data['captured_source_assignments']==2**53+1).all()
 assert (data['unassessed_source_assignments']==7).all() and (data['captured_clones']==1).all()
 parent=pl.read_parquet(out/'S1/gel.parquet');assert parent.filter(pl.col('clone_row')==1)['any_variant_covered'].null_count()==3
 print('PASS: parent counts exact beyond 2^53, no variant duplication, missing capture stays null')
