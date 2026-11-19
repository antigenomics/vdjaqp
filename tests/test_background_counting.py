"""Excluded references must not remove a residue they never contributed."""
import numpy as np
import polars as pl
import pytest
from vdjaqp.background import Background


def test_background_subtracts_only_counted_slots():
    background = Background.__new__(Background)
    background.code = np.array([0] * 38 + [18])
    background.length = np.array([38, 1])
    background.cell_id = np.zeros(39, dtype=int)
    background.n_cells = 1
    backoff = np.full(19, .999 / 18)
    backoff[18] = .001
    background._backoff = lambda counts, pooled: backoff[None, :]
    background._tabulate(pl.DataFrame({'clone_id': ['candidate', 'Ab_reference']}), 'Ab_')
    assert background._contributed.tolist() == [True] * 38 + [False]
    letters = np.arange(19)
    included, excluded = np.zeros(19, dtype=int), np.full(19, 38)
    substitution = (np.ones((19, 19)) - np.eye(19)) / 18
    frequency = background.frequency(excluded, letters)
    assert np.isclose(frequency.sum(), 1, rtol=0, atol=3e-16)
    assert np.allclose(background.mixed(excluded, letters, substitution), frequency @ substitution, rtol=1e-15, atol=1e-17)
    old_frequency = np.maximum(background.numerator[0] - (letters == 0), 1e-9) / (background.total[0] - 1.0 + 20.0)
    old_mixed = np.maximum((background.numerator @ substitution)[0] - substitution[0], 1e-9) / (background.total[0] - 1.0 + 20.0)
    assert np.array_equal(background.frequency(included, letters), old_frequency)
    assert np.array_equal(background.mixed(included, letters, substitution), old_mixed)


def test_explicit_membership_survives_full_population_rename(monkeypatch):
    from vdjaqp import background as bg
    from vdjaqp.airr_ms_inputs import attach_background_membership
    from vdjaqp.evidence import AlignmentScorer
    from vdjaqp.errormodel import ErrorModel,mass_prior
    from vdjaqp.aggregate import ClonotypeScorer
    n=33052
    ids=[f'Ab_ref{i}' if i<12 else f'candidate{i}' for i in range(n)]
    frame=pl.DataFrame(dict(clone_id=ids,seq=['ACY']*12+['ACD']*(n-12),
        region_ids=[[6]*3]*n,cdr3_start=[0]*n,cdr3_end=[3]*n,
        v_name=['IGHV1']*n,j_name=['IGHJ1']*n,c_name=['IGHG1']*n))
    def anchor(self,clones):
        self.cells=pl.DataFrame(dict(cell_id=[0],side=['N'],gene=[''],germline_pos=[-1]))
        return np.zeros(self.code.size,dtype=np.int32),1
    monkeypatch.setattr(bg.Background,'_anchor',anchor)
    tagged=attach_background_membership(frame,inserted_reference_ids=ids[:12])
    assert tagged.drop('background_member').equals(frame)
    renamed=tagged.with_columns(pl.Series('clone_id',[f'neutral{i}' for i in range(n)]))
    legacy=Background(frame);explicit=Background(tagged);neutral=Background(renamed)
    for other in [explicit,neutral]:
        for name in ['total','pooled','numerator','_contributed','cell_id']:
            np.testing.assert_array_equal(getattr(legacy,name),getattr(other,name))
        slots=np.repeat(np.array([0,36]),19);letters=np.tile(np.arange(19),2)
        np.testing.assert_array_equal(legacy.frequency(slots,letters),other.frequency(slots,letters))
        np.testing.assert_array_equal(legacy.mixed(slots,letters,mass_prior()),other.mixed(slots,letters,mass_prior()))
    assert neutral.total.sum()==(n-12)*3
    # Demonstrate the erased-role path actually changes the background.
    wrong=Background(renamed.drop('background_member'))
    assert wrong.total.sum()==n*3 and not np.array_equal(wrong.pooled,neutral.pooled)
    peptides=pl.DataFrame(dict(peptide_id=['p'],spectrum_id=['s'],seq=['ACY'],
        local_confidence=[[.9]*3],search_type=['db'],alc=[.9],alc_rank=[1.],n_reads=[1]))
    hits=pl.DataFrame(dict(peptide_id=['p','p'],clone_id=[ids[0],ids[12]],
                           start=[0,0],end=[3,3],n_subs=[0,1]))
    error=ErrorModel(np.array([0.,1.]),np.array([.1,.01]),mass_prior(),.1)
    a=AlignmentScorer(tagged,error,explicit).score(hits,peptides)
    renamed_hits=hits.with_columns(pl.Series('clone_id',['neutral0','neutral12']))
    b=AlignmentScorer(renamed,error,neutral).score(renamed_hits,peptides)
    np.testing.assert_array_equal(a['llr'].to_numpy(),b['llr'].to_numpy())
    x=ClonotypeScorer(tagged,explicit).score(a,peptides)
    y=ClonotypeScorer(renamed,neutral).score(b,peptides).with_columns(
        pl.col('clone_id').replace_strict(dict(zip(renamed['clone_id'],ids))))
    assert x.sort('clone_id').equals(y.sort('clone_id'))


@pytest.mark.parametrize('values', [[True,None],[1,0],['true','false']])
def test_background_member_requires_boolean_nonmissing(values):
    background=Background.__new__(Background)
    background.code=np.array([0,1]);background.length=np.array([1,1])
    with pytest.raises(ValueError,match='background_member'):
        background._tabulate(pl.DataFrame(dict(clone_id=['a','b'],background_member=values)),'Ab_')
