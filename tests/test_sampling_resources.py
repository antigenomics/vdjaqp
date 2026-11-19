"""Invented counts exercise worker sizing without stochastic/scientific draws."""
import importlib.util
from pathlib import Path
import sys
import numpy as np
import polars as pl
import pytest


@pytest.mark.parametrize('requested,budget,expected', [('all',2,[2]),('9',2,[2]),('1',1,[])])
def test_sampling_pool_uses_resolved_allocation(tmp_path, monkeypatch, requested,budget,expected):
    import vdjaqp
    path=Path(__file__).parents[1]/'scripts/simulate_repertoire_sampling.py'
    spec=importlib.util.spec_from_file_location('sampling_resources',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    source=tmp_path/'abundance.parquet';processing=tmp_path/'processing.tsv'
    pl.DataFrame({'corrected_junction_nt_key':['a','b']*3,
        'subset_id':['L1','L1','L2','L2','L3','L3'],
        'molecule_groups':[2,1,2,1,2,1]}).write_parquet(source)
    pl.DataFrame({'input':['L1','L2','L3'],'molecule_groups':[4,4,4]}).write_csv(processing,separator='\t')
    budgets=[];workers=[]
    monkeypatch.setattr(vdjaqp,'resolve_cpu_budget',lambda request: budgets.append(request) or budget)
    def fake_draw(task):
        counts, depth, replicates, seed = task
        return np.tile(counts, (replicates, 1)).astype(np.uint32)
    monkeypatch.setattr(module, 'draw_library', fake_draw)
    class Pool:
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def map(self,fn,tasks):
            return [np.tile(counts,(replicates,1)).astype(np.uint32)
                for counts,depth,replicates,seed in tasks]
    class Context:
        def Pool(self,size):workers.append(size);return Pool()
    monkeypatch.setattr(module.mp,'get_context',lambda method: Context())
    monkeypatch.setattr(module.mp,'cpu_count',lambda: 64)
    monkeypatch.setattr(sys,'argv',['sampling','--abundance',str(source),
        '--processing',str(processing),'--output',str(tmp_path/'out'),'--cpus',requested])
    module.main()
    assert budgets==[requested] and workers==expected
