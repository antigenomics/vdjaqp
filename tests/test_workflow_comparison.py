"""A partial or mismatched table cannot masquerade as full reproduction."""
from pathlib import Path
import subprocess
import sys
import polars as pl


def test_rejects_missing_endpoint_and_wrong_context(tmp_path):
    columns='n_candidates n_positives n_benchmark_negatives average_precision roc_auc threshold_confirmed_nats true_positive false_positive false_negative true_negative precision_at20 recall_at20 f1_at20 n_clonotypes median_rank worst_rank top_50 top_500'.split()
    rows=[]
    for case in ['gel','solution',*[f'individual_{i:02d}' for i in range(1,9)]]:
        for method in ['alignment_union','exclusive_spectrum']:
            row=dict.fromkeys(columns,0)
            row.update(case=case,method=method,n_candidates=33052,n_clonotypes=33052,
                n_positives=12 if case in ['gel','solution'] else 1,
                threshold_confirmed_nats=20,average_precision=.5,roc_auc=.5)
            rows.append(row)
    fresh=pl.DataFrame(rows)
    accepted=pl.concat([fresh.filter(pl.col('method')=='alignment_union').with_columns(method=pl.lit(label)) for label in ['original arithmetic alignment_union','corrected arithmetic alignment_union']]+[fresh.filter(pl.col('method')=='exclusive_spectrum').with_columns(method=pl.lit('corrected arithmetic exclusive_spectrum'))])
    old=tmp_path/'accepted.tsv';accepted.write_csv(old,separator='\t')
    script=Path(__file__).resolve().parents[1]/'scripts/compare_workflow_reproduction.py'
    invalid=[fresh.drop('average_precision'),fresh.with_columns(pl.col('case').replace({'individual_08':'wrong_context'})),fresh.with_columns(roc_auc=pl.lit(float('nan')))]
    for i,frame in enumerate(invalid):
        source=tmp_path/f'bad{i}.tsv';frame.write_csv(source,separator='\t')
        result=subprocess.run([sys.executable,str(script),str(source),str(old),str(tmp_path/f'out{i}')],capture_output=True,text=True,timeout=10)
        assert result.returncode!=0
        assert not (tmp_path/f'out{i}').exists()
