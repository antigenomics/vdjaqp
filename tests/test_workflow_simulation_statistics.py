"""Post-scoring reference metrics use exported canonical IDs, not source labels."""
from pathlib import Path
import runpy
import sys
from types import SimpleNamespace

import polars as pl


def test_canonical_panel_ids_after_lossless_private_mapping(tmp_path,monkeypatch):
    scripts=Path(__file__).parents[1]/'scripts'
    monkeypatch.syspath_prepend(str(scripts))
    import simulate_airr_ms
    import compress_airr_ms_artifacts
    import vdjaqp
    repertoire=tmp_path/'source.parquet';panel=tmp_path/'source_panel.csv';out=tmp_path/'final';out.mkdir()
    source=[f'source_{i}' for i in range(13)];canonical=[f'record_{i}' for i in range(13)]
    pl.DataFrame(dict(clone_id=source,background_member=[False]*12+[True])).write_parquet(repertoire)
    pl.DataFrame(dict(clone_id=source[:12])).write_csv(panel)
    pl.DataFrame(dict(chain_id=canonical)).write_parquet(out/'physical_truth.parquet')
    pl.DataFrame(dict(chain_id=canonical[:12])).write_parquet(out/'pooled_panel_coverage.parquet')
    scores=pl.DataFrame(dict(clone_id=canonical,score=list(range(13)),confirmed_nats=[0.]*13))
    for enzyme in ('Trypsin','Chymotrypsin','pooled'):scores.write_parquet(out/(enzyme+'_scores.parquet'))
    monkeypatch.setattr(simulate_airr_ms,'parse_args',lambda:SimpleNamespace(cpus=2,repertoire=repertoire,target_ids=panel,output=out))
    monkeypatch.setattr(simulate_airr_ms,'main',lambda:None)
    monkeypatch.setattr(compress_airr_ms_artifacts,'compress_case',lambda root:None)
    monkeypatch.setattr(vdjaqp,'resolve_cpu_budget',lambda cpus:cpus)
    runpy.run_path(str(scripts/'workflow_simulation.py'),run_name='__main__')
    statistics=pl.read_csv(out/'classification_statistics.tsv',separator='\t')
    assert statistics.height==6
    assert statistics.filter(pl.col('metric_scope')=='reference_retrieval')['n_positives'].to_list()==[12]*3
    assert statistics.filter(pl.col('metric_scope')=='physical_record_presence')['average_precision'].null_count()==3
