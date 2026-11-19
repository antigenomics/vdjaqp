"""Synthetic regression checks for numerical boundaries; no scientific inputs."""
from types import SimpleNamespace

import numpy as np
import polars as pl
import pytest

from vdjaqp.alphabet import encode, encode_frame
from vdjaqp.background import Background
from vdjaqp.errormodel import isotonic
from vdjaqp.pgen import sharing_probability


def test_isotonic_pools_ties_with_sample_count_weights():
    x = np.array([.1, .1, .2, .3])
    y = np.array([1., 0., 1., 0.])
    for order in [np.arange(4), np.array([1, 0, 3, 2])]:
        knots, fitted = isotonic(x[order], y[order])
        np.testing.assert_array_equal(knots, [.1, .1, .2, .3])
        np.testing.assert_allclose(fitted, [2/3, 2/3, 2/3, 0.])


@pytest.mark.parametrize('x,y', [([], []), ([0.], []), ([[0.]], [[1.]]),
    ([np.nan], [1.]), ([0.], [np.inf])])
def test_isotonic_rejects_invalid_fit_inputs(x, y):
    with pytest.raises(ValueError, match='finite|nonempty|dimensional'):
        isotonic(np.asarray(x), np.asarray(y))


def test_sharing_probability_boundaries():
    with np.errstate(all='raise'):
        np.testing.assert_array_equal(sharing_probability([0., .5, 1.], 0), [0., 0., 0.])
        np.testing.assert_allclose(sharing_probability([0., .5, 1.], 2), [0., .75, 1.])


@pytest.mark.parametrize('probability,depth', [([-1.], 1), ([1.1], 1),
    ([np.nan], 1), ([np.inf], 1), ([.1], -1), ([.1], np.nan), ([.1], np.inf)])
def test_sharing_probability_rejects_invalid_inputs(probability, depth):
    with pytest.raises(ValueError):
        sharing_probability(probability, depth)


@pytest.mark.parametrize('code,membership', [([0], [False]), ([-1], [True])])
def test_background_requires_usable_population(code, membership):
    background = Background.__new__(Background)
    background.code = np.array(code)
    background.length = np.array([1])
    background.cell_id = np.array([0])
    background.n_cells = 1
    with pytest.raises(ValueError, match='background.*canonical'):
        background._tabulate(pl.DataFrame({'clone_id': ['c'], 'background_member': membership}), 'Ab_')


def test_encode_frame_unicode_offsets_match_single_sequence_encoding():
    sequences = ['AéC', 'D😀E', 'FG', '']
    codes, offsets, lengths = encode_frame(pl.Series(sequences))
    np.testing.assert_array_equal(lengths, [3, 3, 2, 0])
    for seq, offset, length in zip(sequences, offsets, lengths):
        np.testing.assert_array_equal(codes[offset:offset+length], encode(seq))


def test_clonotype_table_preserves_undefined_ndn_coverage(monkeypatch):
    from vdjaqp import report
    clones = pl.DataFrame({'clone_id': ['c'], 'seq': ['ACFW'], 'cdr3_seq': ['CFW'],
        'v_name': ['IGHV1'], 'j_name': ['IGHJ1'], 'c_name': ['IGHG1'],
        'lineage_id': ['l'], 'nt_seq': ['GCTTGTTTTTGG']})
    profile = SimpleNamespace(regions=np.array([0, 5, 5, 7]), seq='ACFW',
        shm=np.zeros(4, dtype=bool), peptides=pl.DataFrame(), is_probe=False,
        covered=lambda **kwargs: float('nan'))
    monkeypatch.setattr(report, 'profile', lambda *args, **kwargs: profile)
    scores = pl.DataFrame({'clone_id': ['c'], 'score': [0.], 'confirmed_nats': [0.],
        'junction_nats': [0.], 'framework_only': [False]})
    result = report.clonotype_table(SimpleNamespace(clonotypes=clones), clone_ids=["c"], scored=scores, targets=set())
    value = result['ndn_covered_pct'][0]
    assert value is None or np.isnan(value)
