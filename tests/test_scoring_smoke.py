"""Bounded synthetic scoring and model-API regression checks; no cohort runs."""
import numpy as np
import polars as pl
import pytest
from vdjaqp.pgen import sharing_probability, junction_pgen


def test_sharing_probability_is_monotone_in_depth():
    """Retained data-independent legacy regression."""
    p = np.array([1e-8])
    assert sharing_probability(p, 1e3) < sharing_probability(p, 1e6)


def test_bundled_pgen_batch_schema_and_anchor_convention():
    frame = junction_pgen(['CARGGDYW', 'ARGGDY', 'CARXXW'], threads=1)
    assert frame['junction'].to_list() == ['CARGGDYW', 'ARGGDY', 'CARXXW']
    assert frame['pgen'][0] > 0
    assert frame['pgen'][1] == frame['pgen'][2] == 0
    assert frame['neglog10_pgen'][1] is None
    assert frame['scoreable'].to_list() == [True, False, False]


def test_empty_peptide_export_preserves_candidates(monkeypatch):
    from vdjaqp import background as bg
    from vdjaqp.errormodel import ErrorModel, mass_prior
    from vdjaqp.evidence import AlignmentScorer
    from vdjaqp.aggregate import ClonotypeScorer
    from vdjaqp.match import PeptideMatcher
    monkeypatch.setattr(bg, 'canonical_v', lambda: {})
    monkeypatch.setattr(bg, 'canonical_j', lambda: {})
    clones = pl.DataFrame({'clone_id': ['clone0'], 'seq': ['ACDEFGHKLMNP'],
        'region_ids': [[6]*12], 'cdr3_start': [0], 'cdr3_end': [12],
        'v_name': ['IGHV1-2'], 'j_name': ['IGHJ4'], 'c_name': ['IGHG1']})
    peptides = pl.DataFrame(schema={'peptide_id': pl.String, 'seq': pl.String,
        'local_confidence': pl.List(pl.Float32), 'search_type': pl.String,
        'alc': pl.Float64, 'alc_rank': pl.Float64, 'n_reads': pl.Int64})
    error = ErrorModel(np.array([0., 1.]), np.array([0.1, 0.01]), mass_prior(), 0)
    background = bg.Background(clones)
    hits = PeptideMatcher.from_clonotypes(clones).search([], [], threads=1)
    scored = AlignmentScorer(clones, error, background).score(hits, peptides)
    assert scored.height == 0 and scored.schema['llr'] == pl.Float64
    result = ClonotypeScorer(clones, background).score(scored, peptides)
    assert result['clone_id'].to_list() == ['clone0']
    assert result['confirmed_nats'].to_list() == [0.0]
    assert np.isfinite(result['score'].to_numpy()).all()


@pytest.mark.parametrize("isotypes", [["IGHG1", "IGHG1"], ["IGHM", "IGHM"], ["IGHG1", "IGHM"]])
@pytest.mark.parametrize("missing_db_confidence", [False, True])
def test_search_to_scoring_with_explicit_error_model(monkeypatch, isotypes, missing_db_confidence):
    from vdjaqp import background as bg
    from vdjaqp.errormodel import ErrorModel, mass_prior
    from vdjaqp.evidence import AlignmentScorer
    from vdjaqp.aggregate import ClonotypeScorer
    from vdjaqp.match import PeptideMatcher
    # Explicit synthetic coordinates avoid reference acquisition and calibration.
    monkeypatch.setattr(bg, 'canonical_v', lambda: {})
    monkeypatch.setattr(bg, 'canonical_j', lambda: {})
    clones = pl.DataFrame({'clone_id': ['clone0', 'clone1'],
        'seq': ['ACDEFGHKLMNP', 'QRSTVWYACDEF'],
        'region_ids': [[6]*12, [6]*12], 'cdr3_start': [0,0], 'cdr3_end': [12,12],
        'v_name': ['IGHV1-2']*2, 'j_name': ['IGHJ4']*2, 'c_name': isotypes})
    peptides = pl.DataFrame({'peptide_id': ['run/1000/2/500/10/1/ACDEFG'],
        'seq': ['ACDEFG'], 'local_confidence': [[float('nan') if missing_db_confidence else 0.99]*6], 'search_type': ['db'],
        'alc': [float('nan')], 'alc_rank': [1.0], 'n_reads': [1]})
    error = ErrorModel(np.array([0.,1.]), np.array([0.1,0.01]), mass_prior(), 0)
    background = bg.Background(clones)
    hits = PeptideMatcher.from_clonotypes(clones).search(peptides['peptide_id'], peptides['seq'], max_subs=0, threads=1)
    scored = AlignmentScorer(clones, error, background).score(hits, peptides)
    result = ClonotypeScorer(clones, background).score(scored, peptides)
    assert set(result['clone_id']) == {'clone0', 'clone1'}
    assert np.isfinite(result['score'].to_numpy()).all()
    assert result['clone_id'][0] == 'clone0'
    assert result['junction_nats'][0] > 0
    assert not result['framework_only'].any()


@pytest.mark.parametrize('chunk_size', [1, 2, 3, 100])
@pytest.mark.parametrize('same_spectrum, expected', [(True, .8), (False, .64)])
def test_exclusive_spectrum_coverage(monkeypatch, chunk_size, same_spectrum, expected):
    """Categorical alternatives sum; distinct observations form a union."""
    from types import SimpleNamespace
    from vdjaqp import aggregate
    scorer = object.__new__(aggregate.ClonotypeScorer)
    scorer.background = SimpleNamespace(code=np.array([0, 1, 0, 1]),
                                        offset=np.array([0, 2]))
    scorer._clone_index = pl.DataFrame({'clone_id': ['a', 'b']})
    # Interleave spectra and include a second candidate and mismatching residue.
    keyed = pl.DataFrame({'clone_row': [0, 1, 0], 'llr': [2., -2., 2.],
        'pep_row': [0, 1, 0], 'start': [0, 0, 0],
        'spectrum_id': ['z', 'a', 'z' if same_spectrum else 'b']})
    monkeypatch.setattr(aggregate, '_CHUNK', chunk_size)
    args = (keyed, np.array([.4, .25, .4]), np.array([0, 1, 0, 2]),
            np.array([0, 2]), np.array([2, 2]))
    confirmed, contradiction = scorer._cover(*args, aggregation='exclusive_spectrum')
    np.testing.assert_allclose(confirmed, [expected, expected, .25, 0])
    np.testing.assert_allclose(contradiction, [0, -.5])
    old, old_contradiction = scorer._cover(*args)
    np.testing.assert_allclose(old, [.64, .64, .25, 0])
    np.testing.assert_array_equal(contradiction, old_contradiction)


def test_aggregation_api_single_alternative_no_competition_and_empty(monkeypatch):
    """Public option preserves deterministic ablation and default single-hit case."""
    from types import SimpleNamespace
    from vdjaqp import aggregate
    clones = pl.DataFrame({'clone_id': ['a'], 'seq': ['AC'],
        'region_ids': [[6, 6]], 'v_name': ['IGHV1-2'], 'j_name': ['IGHJ4'],
        'c_name': ['IGHG1']})
    background = SimpleNamespace(code=np.array([0, 1]), offset=np.array([0]),
                                 information=lambda: np.ones(2))
    scorer = aggregate.ClonotypeScorer(clones, background)
    peptides = pl.DataFrame({'peptide_id': ['p', 'q'], 'seq': ['AC', 'AC'],
        'spectrum_id': ['s', 's'], 'n_reads': [1, 1], 'alc': [.5, .5],
        'alc_rank': [1., 1.]})
    scored = pl.DataFrame({'clone_id': ['a', 'a'], 'peptide_id': ['p', 'q'],
                          'start': [0, 0], 'llr': [2., 2.]})
    for alignments, competition in [(scored.head(1), True), (scored, False),
                                    (scored.head(0), True)]:
        old = scorer.score(alignments, peptides, competition=competition)
        new = scorer.score(alignments, peptides, competition=competition,
                           aggregation='exclusive_spectrum')
        assert old.equals(new)
    old = scorer.score(scored, peptides)
    new = scorer.score(scored, peptides, aggregation='exclusive_spectrum')
    assert new['confirmed_nats'].item() > old['confirmed_nats'].item()
    with pytest.raises(ValueError, match='aggregation'):
        scorer.score(scored, peptides, aggregation='unknown')


@pytest.mark.parametrize('chunk_size', [1, 7, 100])
def test_exclusive_coverage_matches_categorical_reference(monkeypatch, chunk_size):
    """Mixed candidates, offsets and wholly disagreeing spectra use the same law."""
    from types import SimpleNamespace
    from vdjaqp import aggregate
    scorer = object.__new__(aggregate.ClonotypeScorer)
    scorer.background = SimpleNamespace(code=np.array([0, 1, 0, 1]),
                                        offset=np.array([0, 2]))
    scorer._clone_index = pl.DataFrame({'clone_id': ['a', 'b']})
    index = np.arange(31)
    keyed = pl.DataFrame({'clone_row': index % 2, 'llr': np.zeros(31),
        'pep_row': index % 3, 'start': np.zeros(31, dtype=int),
        'spectrum_id': [str(i % 5) for i in index]})
    gamma = (1 + index % 3) / 100
    pep_code = np.array([0, 1, 0, 2, -1, -1])
    # Small independent scalar oracle: no residue expansion/grouping helpers.
    mass = {}
    for row, g in zip(keyed.iter_rows(named=True), gamma):
        for residue in range(2):
            aa = pep_code[2 * row['pep_row'] + residue]
            if aa >= 0 and aa == scorer.background.code[2 * row['clone_row'] + residue]:
                key = (row['spectrum_id'], 2 * row['clone_row'] + residue)
                mass[key] = mass.get(key, 0.) + g
    miss = np.ones(4)
    for (_, slot), value in mass.items():
        miss[slot] *= 1 - value
    monkeypatch.setattr(aggregate, '_CHUNK', chunk_size)
    observed, _ = scorer._cover(keyed, gamma, pep_code, np.array([0, 2, 4]),
        np.array([2, 2, 2]), aggregation='exclusive_spectrum')
    np.testing.assert_allclose(observed, 1 - miss, rtol=1e-14, atol=1e-15)
    # No agreeing residue in any spectrum must also handle empty sparse groups.
    empty, _ = scorer._cover(keyed, gamma, np.full(6, -1), np.array([0, 2, 4]),
        np.array([2, 2, 2]), aggregation='exclusive_spectrum')
    np.testing.assert_array_equal(empty, np.zeros(4))
