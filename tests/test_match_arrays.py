"""Flat exports preserve native hit order, coordinates and buffer ownership."""
import gc

import polars as pl
import pytest

from vdjaqp.match import PeptideMatcher


def test_sequence_aliases_preserve_native_outputs():
    matcher = PeptideMatcher(['a', 'b'], ['ACDEFGHIKLMNPQACDEFGHIK', 'ACDEYGHIKLMNPQ'])
    ids = ['same', 'other', 'same', 'il', 'short', 'miss']
    seqs = ['ACDEFGHIK', 'ACDEYGHIK', 'ACDEFGHIK', 'ACDEFGHLK', 'ACD', 'WWWWWWWWW']
    for indels in (0, 1):
        actual = matcher.search(ids, seqs, max_subs=0 if indels else 1, max_indels=indels, threads=1)
        # One query at a time is an independent witness of expansion and ordering.
        import polars as pl
        expected = pl.concat([matcher.search([i], [s], max_subs=0 if indels else 1, max_indels=indels, threads=1)
                              for i, s in zip(ids, seqs)])
        assert actual.equals(expected)


def test_budget_groups_preserve_alias_order():
    import polars as pl
    matcher = PeptideMatcher(['a'], ['ACDEFGHIKLMNPQRSTVWYACDEFGHIK'])
    ids = ['long', 'short', 'alias', 'long2']
    seqs = ['ACDEFGHIKLMNPQ', 'ACDEFGHIK', 'ACDEFGHLK', 'ACDEFGHIKLMNPQ']
    budget = lambda n: n // 6
    actual = matcher.search(ids, seqs, budget=budget, threads=1)
    order = sorted(range(len(ids)), key=lambda i: budget(len(seqs[i])))
    expected = pl.concat([matcher.search([ids[i]], [seqs[i]], budget=budget, threads=1) for i in order])
    assert actual.equals(expected)


@pytest.mark.parametrize('indels', [0, 1])
def test_flat_export_matches_native_hit_objects(indels):
    """The independent oracle reads native objects, never the flat export."""
    matcher = PeptideMatcher(
        ['c0', 'c1'], ['ACDEFGHIKLMNPQACDEFGHIK', 'ACYEFGHIKLMNPQ'])
    ids = ['long', 'short', 'repeated', 'il', 'absent', 'too_short']
    seqs = ['ACDEFGHIKLMNPQ', 'ACDEFGHIK', 'ACDEFGHIK',
            'ACDEFGHLK', 'WWWWWWWWW', 'ACD']
    budget_for = lambda n: n // 6 - indels
    actual = matcher.search(ids, seqs, budget=budget_for,
                            max_indels=indels, threads=1)
    rows = []
    for budget in sorted({budget_for(len(s)) for s in seqs if len(s) >= 6}):
        indices = [i for i, s in enumerate(seqs)
                   if len(s) >= 6 and budget_for(len(s)) == budget]
        native = matcher._index.search_batch(
            [seqs[i].replace('I', 'L') for i in indices],
            max_subs=budget, max_indels=indels, threads=1)
        for query, i in enumerate(indices):
            for hit in native[query]:
                rows.append((matcher.clone_ids[hit.ref_id], ids[i], hit.offset,
                             hit.offset + len(seqs[i]), hit.n_subs,
                             [int(pos) for pos, _, _ in hit.mismatches]))
    expected = pl.DataFrame(rows, schema=actual.schema, orient='row')
    del native, matcher
    gc.collect()
    assert actual.equals(expected)
    assert actual.filter(pl.col('peptide_id') == 'short').height > 1
    if not indels:
        assert [2] in actual['mismatch_positions'].to_list()


def test_empty_flat_export_keeps_schema():
    matcher = PeptideMatcher(['c0'], ['ACDEFGHIKLMNPQ'])
    empty = matcher.search([], [], threads=1)
    no_hits = matcher.search(['q0'], ['WWWWWWWWW'], max_subs=0, threads=1)
    assert empty.is_empty() and no_hits.is_empty()
    assert empty.schema == no_hits.schema
    assert empty.schema['mismatch_positions'] == pl.List(pl.Int16)


@pytest.mark.parametrize('field,value', [('offset', 2**31), ('n_subs', 2**15)])
def test_flat_export_rejects_out_of_schema_values(field, value):
    """Signed export limits must raise, as the original strict table did."""
    import numpy as np

    class Result:
        def to_numpy(self):
            arrays = {
                'query_begin': np.array([0, 1], dtype=np.uint32),
                'ref_id': np.array([0], dtype=np.uint32),
                'offset': np.array([0], dtype=np.uint32),
                'n_subs': np.array([0], dtype=np.uint16),
                'mm_begin': np.array([0, 0], dtype=np.uint32),
                'mm_pos': np.array([], dtype=np.uint16),
            }
            arrays[field][0] = value
            return arrays

    class Index:
        def search_batch(self, *args, **kwargs):
            return Result()

    matcher = PeptideMatcher.__new__(PeptideMatcher)
    matcher.clone_ids = ['c0']
    matcher._index = Index()
    with pytest.raises(pl.exceptions.InvalidOperationError):
        matcher.search(['q0'], ['ACDEFGHIK'], threads=1)


def test_flat_export_preserves_separate_observation_scoring(monkeypatch):
    import numpy as np
    from vdjaqp import background as bg
    from vdjaqp.aggregate import ClonotypeScorer
    from vdjaqp.errormodel import ErrorModel, mass_prior
    from vdjaqp.evidence import AlignmentScorer

    monkeypatch.setattr(bg, 'canonical_v', lambda: {})
    monkeypatch.setattr(bg, 'canonical_j', lambda: {})
    clones = pl.DataFrame({
        'clone_id': ['c0', 'c1'],
        'seq': ['ACDEFGHIKLMNPQ', 'ACYEFGHIKLMNPQ'],
        'region_ids': [[6]*14]*2, 'cdr3_start': [0, 0], 'cdr3_end': [14, 14],
        'v_name': ['IGHV1-2']*2, 'j_name': ['IGHJ4']*2,
        'c_name': ['IGHG1', 'IGHM']})
    peptides = pl.DataFrame({
        'peptide_id': ['p0', 'p1'], 'spectrum_id': ['s0', 's1'],
        'seq': ['ACDEFGHIK', 'ACDEFGHIK'],
        'local_confidence': [[.99]*9, [.65]*9], 'search_type': ['denovo']*2,
        'alc': [.99, .65], 'alc_rank': [1., 2.], 'n_reads': [1, 2]})
    actual = PeptideMatcher.from_clonotypes(clones).search(
        peptides['peptide_id'], peptides['seq'], max_subs=1, threads=1)
    expected = pl.DataFrame([
        ('c0', 'p0', 0, 9, 0, []), ('c1', 'p0', 0, 9, 1, [2]),
        ('c0', 'p1', 0, 9, 0, []), ('c1', 'p1', 0, 9, 1, [2]),
    ], schema=actual.schema, orient='row')
    assert actual.equals(expected)
    background = bg.Background(clones)
    error = ErrorModel(np.array([0., 1.]), np.array([.1, .01]), mass_prior(), .1)
    scorer = AlignmentScorer(clones, error, background)
    scored = scorer.score(actual, peptides)
    reference = scorer.score(expected, peptides)
    assert scored.equals(reference)
    assert scored['llr'][0] != scored['llr'][2]
    aggregate = ClonotypeScorer(clones, background)
    for mode in ['alignment_union', 'exclusive_spectrum']:
        assert aggregate.score(scored, peptides, aggregation=mode).equals(
            aggregate.score(reference, peptides, aggregation=mode))
