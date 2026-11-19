"""Unsafe gapped geometry and unpaired identifiers must fail before scoring."""
import pytest
from vdjaqp.match import PeptideMatcher


@pytest.mark.parametrize('query,counts', [('ACDDEFGHIK',(1,0)), ('ACDEGHIK',(0,1))])
def test_gapped_search_exports_native_reference_span(query, counts):
    matcher = PeptideMatcher(['c'], ['ACDEFGHIK'])
    matches = matcher.search(['p'], [query], max_subs=0, max_indels=1, threads=1)
    assert matches.select('start','end','n_ins','n_dels').rows()==[(0,9,*counts)]
    from vdjaqp.evidence import AlignmentScorer
    import polars as pl
    scorer=AlignmentScorer.__new__(AlignmentScorer)
    with pytest.raises(ValueError, match='gapped'):
        scorer.score(matches, pl.DataFrame())


def test_reference_identifiers_must_match_sequences():
    with pytest.raises(ValueError, match='same length|aligned'):
        PeptideMatcher(['c', 'missing'], ['ACDEFGH'])


def test_query_identifiers_must_match_sequences_even_without_usable_queries():
    matcher = PeptideMatcher(['c'], ['ACDEFGH'])
    with pytest.raises(ValueError, match='same length|aligned'):
        matcher.search([], ['A'], threads=1)
