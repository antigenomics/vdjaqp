"""A valid native evidence table need not carry archival entropy annotations."""
from types import SimpleNamespace
import numpy as np
import polars as pl
from vdjaqp.report import profile


def test_profile_missing_entropy_is_unavailable_without_losing_placements(monkeypatch):
    monkeypatch.setattr('vdjaqp.report.shm_positions', lambda seq, *args: np.zeros(len(seq), dtype=bool))
    c = pl.DataFrame(dict(clone_id=['c'], seq=['CAW'], region_ids=[[4,6,7]],
                         cdr3_start=[0], cdr3_end=[3], v_name=['V'], j_name=['J'], c_name=['unknown']))
    ev = SimpleNamespace(clonotypes=c,
        alignments=pl.DataFrame(dict(clone_id=['c'], peptide_id=['p'], start=[0], end=[3], llr=[1.])),
        peptides=pl.DataFrame(dict(peptide_id=['p'],seq=['CAW'],alc=[90.],search_type=['x'],enzyme=['x'])),
        background=SimpleNamespace(offset=np.array([0]),information=lambda:np.ones(3)))
    result = profile('c', ev, targets=set())
    assert np.isnan(result.entropy).all() and result.entropy.shape == (3,)
    assert result.agreeing_depth.tolist() == [1,1,1]
    assert result.peptides.height == 1


def test_explicit_nucleotide_support_preserves_observed_constant_residues(monkeypatch):
    monkeypatch.setattr('vdjaqp.report.shm_positions', lambda seq, *args: np.zeros(len(seq), dtype=bool))
    c = pl.DataFrame(dict(clone_id=['c'], seq=['CAWG'], region_ids=[[4,6,7,8]],
                         cdr3_start=[0], cdr3_end=[4], v_name=['V'], j_name=['J'], c_name=['unknown']))
    ev = SimpleNamespace(clonotypes=c,
        alignments=pl.DataFrame(dict(clone_id=['c'], peptide_id=['p'], start=[0], end=[4], llr=[1.])),
        peptides=pl.DataFrame(dict(peptide_id=['p'],seq=['CAWG'],alc=[90.],search_type=['x'],enzyme=['x'])),
        background=SimpleNamespace(offset=np.array([0]),information=lambda:np.ones(4)))
    assert profile('c', ev, targets=set()).nt_support.tolist() == [True,True,True,False]
    ev.clonotypes=c.with_columns(pl.Series('nt_support',[[True,True,True,True]]))
    assert profile('c', ev, targets=set()).nt_support.tolist() == [True,True,True,True]
    ev.clonotypes=c.with_columns(pl.Series('nt_support',[[True]]))
    import pytest
    with pytest.raises(ValueError,match='nucleotide support'):
        profile('c', ev, targets=set())
