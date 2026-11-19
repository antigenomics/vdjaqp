"""Equivalent gene aliases must retain their shared measured-sequence context."""
import numpy as np
import polars as pl
import pytest

from vdjaqp import background, shm
from vdjaqp.germline import common_reference


@pytest.mark.parametrize('call', ['V1,V2', ' V2,V1,V2 ', 'V1*01,V2*02'])
def test_common_reference_requires_every_alternative(call):
    assert common_reference(call, {'V1': 'ACDE', 'V2': 'ACDE'}) == ('V1,V2', 'ACDE')
    assert common_reference(call, {'V1': 'ACDE'}) is None
    assert common_reference(call, {'V1': 'ACDE', 'V2': 'ACDF'}) is None


def test_equivalent_aliases_match_single_gene_background_and_shm(monkeypatch):
    for module in (background, shm):
        monkeypatch.setattr(module, 'canonical_v', lambda: {'V1-1': 'ACDE', 'V1-2': 'ACDE'})
        monkeypatch.setattr(module, 'canonical_j', lambda: {'J1': 'FGHW', 'J2': 'FGHW'})
    frame = pl.DataFrame(dict(clone_id=['x'], seq=['ACDEFGHW'],
        region_ids=[[0, 0, 5, 5, 6, 6, 7, 7]], cdr3_start=[2], cdr3_end=[6],
        v_name=['V1-1'], j_name=['J1'], c_name=['C']))
    alternatives = frame.with_columns(v_name=pl.lit('V1-2,V1-1,V1-2'), j_name=pl.lit('J2,J1'))
    single, aliases = background.Background(frame), background.Background(alternatives)
    assert single.cells.drop('gene').equals(aliases.cells.drop('gene'))
    for name in ('cell_id', 'numerator', 'total', 'pooled'):
        np.testing.assert_array_equal(getattr(single, name), getattr(aliases, name))
    np.testing.assert_array_equal(single.information(), aliases.information())
    assert shm.revert_junctions(frame).equals(shm.revert_junctions(alternatives))
    assert alternatives['v_name'][0] == 'V1-2,V1-1,V1-2'


def test_shared_reference_does_not_assign_ambiguous_counts_to_first_gene(monkeypatch):
    monkeypatch.setattr(background,'canonical_v',lambda: {'IGHV1-1':'ACD','IGHV1-2':'ACD'})
    monkeypatch.setattr(background,'canonical_j',lambda: {})
    frame=pl.DataFrame(dict(clone_id=['a','b','ambiguous'],seq=['ACD','CCD','DCD'],
        region_ids=[[0]*3]*3,cdr3_start=[3]*3,cdr3_end=[3]*3,
        v_name=['IGHV1-1','IGHV1-2','IGHV1-2,IGHV1-1'],j_name=['missing']*3,c_name=['C']*3))
    model=background.Background(frame)
    first=model.cells.filter((pl.col('side')=='V') & (pl.col('germline_pos')==0))
    assert set(first['gene']) == {'IGHV1-1','IGHV1-2','IGHV1-1,IGHV1-2'}
    assert model.total[first['cell_id'].to_numpy()].tolist() == [1.,1.,1.]


def test_alternatives_from_different_v_families_use_pooled_backoff():
    model=background.Background.__new__(background.Background)
    model.cells=pl.DataFrame(dict(cell_id=[0,1],side=['V','V'],
        gene=['IGHV1-1,IGHV3-1','IGHV1-1'],germline_pos=[0,0]))
    counts=np.zeros((2,19));counts[0,0]=4;counts[1,1]=6
    pooled=np.full(19,1/19)
    backoff=model._backoff(counts,pooled)
    np.testing.assert_array_equal(backoff[0],(counts[0]+20*pooled)/24)


@pytest.mark.parametrize('missing_side', ['v', 'j'])
def test_one_sided_unresolved_context_preserves_known_side(monkeypatch, missing_side):
    monkeypatch.setattr(background, 'canonical_v', lambda: {'V': 'AC'})
    monkeypatch.setattr(background, 'canonical_j', lambda: {'J': 'DE'})
    # A real aligned position and an insertion on each side; explicit NDN wins.
    monkeypatch.setattr(background, 'aligned_columns',
        lambda query, ref: iter([(0, 0, ref[0])]))
    frame = pl.DataFrame(dict(clone_id=['x'], seq=['ACDEFG'],
        region_ids=[[0, 0, 6, 6, 7, 7]], cdr3_start=[2], cdr3_end=[10],
        v_name=['missing' if missing_side == 'v' else 'V'],
        j_name=['missing' if missing_side == 'j' else 'J'], c_name=['C']))
    model = background.Background(frame)
    if missing_side == 'v':
        assert model.context_valid.tolist() == [False, False, True, True, True, True]
        assert model.cell_id[5] == 0  # Known J insertion remains N.
    else:
        assert model.context_valid.tolist() == [True, True, True, True, False, False]
        assert model.cell_id[1] == 0  # Known V insertion remains N.


@pytest.mark.parametrize('call,expected', [('V2,V1', 'C'), ('V1,missing', 'A'),
                                         ('V1,V3', 'A')])
def test_shm_alternatives_revert_only_resolved_mutated_flank(monkeypatch, call, expected):
    monkeypatch.setattr(shm, 'canonical_v', lambda: {'V1': 'C', 'V2': 'C', 'V3': 'D'})
    monkeypatch.setattr(shm, 'canonical_j', lambda: {'J': 'E'})
    monkeypatch.setattr(shm, 'aligned_columns',
        lambda query, ref: iter([(0 if ref != 'E' else 1, 0, ref)]))
    frame = pl.DataFrame(dict(clone_id=['x'], seq=['AA'], region_ids=[[5, 5]],
        cdr3_start=[0], cdr3_end=[2], v_name=[call], j_name=['J'], c_name=['C']))
    result = shm.revert_junctions(frame)
    assert result['junction_reverted'][0] == expected + 'E'
    assert result['n_reverted'][0] == (2 if expected == 'C' else 1)


def test_shm_display_resolves_identical_alternatives(monkeypatch):
    from vdjaqp import report
    monkeypatch.setattr(report,'canonical_v',lambda: {'V1':'C','V2':'C','V3':'D'})
    monkeypatch.setattr(report,'canonical_j',lambda: {'J':'E'})
    monkeypatch.setattr(report,'aligned_columns',lambda query,ref: iter([(0 if ref!='E' else 1,0,ref)]))
    args=('AA',np.array([5,5]),0,2)
    assert report.shm_positions(args[0],'V2,V1','J',*args[1:]).tolist()==[True,True]
    assert report.shm_positions(args[0],'V1,missing','J',*args[1:]).tolist()==[False,True]
    assert report.shm_positions(args[0],'V1,V3','J',*args[1:]).tolist()==[False,True]


def test_unresolved_framework_is_not_non_templated_evidence(monkeypatch):
    monkeypatch.setattr(background, 'canonical_v', lambda: {})
    monkeypatch.setattr(background, 'canonical_j', lambda: {})
    frame = pl.DataFrame(dict(clone_id=['x'], seq=['ACDE'],
        region_ids=[[0, 6, 7, 8]], cdr3_start=[1], cdr3_end=[2],
        v_name=['missing'], j_name=['missing'], c_name=['C']))
    model = background.Background(frame)
    assert model.context_valid.tolist() == [False, True, False, True]
    assert model.total[0] == 1  # Explicit NDN only, not missing V/J context.
    assert model.total.sum() == 2
    assert np.count_nonzero(model.pooled) == 4  # All observed AA retained globally.
    assert model.information()[[0, 2]].tolist() == [0., 0.]
    from vdjaqp.evidence import AlignmentScorer
    from vdjaqp.errormodel import ErrorModel, mass_prior
    peptides = pl.DataFrame(dict(peptide_id=['p'], seq=['A'],
        local_confidence=[[.99]], search_type=['db'], alc=[.99],
        alc_rank=[1.], n_reads=[1]))
    hits = pl.DataFrame(dict(peptide_id=['p'], clone_id=['x'],
        start=[0], end=[1], n_subs=[0]))
    error = ErrorModel(np.array([0., 1.]), np.array([.1, .01]), mass_prior(), 0)
    assert AlignmentScorer(frame, error, model).score(hits, peptides)['llr'][0] == 0.
