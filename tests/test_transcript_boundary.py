import pytest
from vdjaqp.transcript_boundary import project_transcript_boundaries


def project(qaln='AACCGG', taln='AACCGG', *, p=4, qstart=1, tstart=1,
            query=None, context=None, mask='all', extra_contexts=(), extra_hits=()):
    query = query or qaln.replace('-', '')
    context = context or taln.replace('-', '')
    hit = dict(sequence_id='q', target='t', qstart=qstart,
               qend=qstart + len(qaln.replace('-', '')) - 1,
               tstart=tstart, tend=tstart + len(taln.replace('-', '')) - 1,
               qlen=len(query), tlen=len(context), qaln=qaln, taln=taln)
    return project_transcript_boundaries(
        [dict(sequence_id='q', sequence=query,
              retained_primary_exact_match_mask=[True]*len(query) if mask == 'all' else mask),
         dict(sequence_id='nohit', sequence='AAA')],
        [dict(target='t', sequence=context, mature_start=p), *extra_contexts],
        [hit, *extra_hits]).to_dicts()


def test_matching_mismatching_and_missing_observations():
    for qaln in ('AACCGG', 'AATCGG'):
        rows = project(qaln, 'AACCGG')
        assert rows[0] == dict(sequence_id='q', mature_start=4, status='resolved',
            alignment_count=1, candidate_count=1, supported_alignment_count=1,
            bracketing_alignment_count=1, nonbracketing_alignment_count=0)
        assert rows[1]['mature_start'] is None and rows[1]['status'] == 'no_alignment'
    # Earlier insertion shifts the observed boundary through native projection.
    assert project('AATCCGG', 'AA-CCGG')[0]['mature_start'] == 5
    # Nucleotide coordinate only: there is no codon rounding or base completion.
    assert project(p=3)[0]['mature_start'] == 3


@pytest.mark.parametrize('qaln,taln,kwargs', [
    ('AAC-GG', 'AACCGG', {}),  # first mature base deleted
    ('AA-CGG', 'AACCGG', {}),  # last leader base deleted
    ('AACTCGG', 'AAC-CGG', {}),  # insertion straddles the boundary
    ('CGG', 'CGG', dict(query='AACCGG', context='AACCGG', qstart=4, tstart=4)),
    ('AAC', 'AAC', dict(query='AACCGG', context='AACCGG')),
])
def test_unobserved_or_nonadjacent_markers(qaln, taln, kwargs):
    row = project(qaln, taln, **kwargs)[0]
    assert row['status'] == 'unsupported_boundary' and row['mature_start'] is None
    assert row['nonbracketing_alignment_count'] == 1


def test_both_observed_marker_bases_require_exact_primary_support():
    for mask in (None, [True, True, False, True, True, True],
                 [True, True, True, False, True, True]):
        row = project(mask=mask)[0]
        assert row['mature_start'] is None and row['status'] == 'unsupported_boundary'
        assert row['bracketing_alignment_count'] == 1
        assert row['supported_alignment_count'] == 0
    # Unrelated unsupported positions do not erase the observed marker evidence.
    assert project(mask=[False, True, True, True, True, False])[0]['mature_start'] == 4


def test_ambiguous_query_markers_and_public_context_rejected():
    for query in ('AANCGG', 'AACNGG'):
        row = project(query, 'AACCGG')[0]
        assert row['mature_start'] is None and row['status'] == 'unsupported_boundary'
        assert row['bracketing_alignment_count'] == 1
        assert row['supported_alignment_count'] == 0
    assert project('NACCGG', 'AACCGG')[0]['mature_start'] == 4
    with pytest.raises(ValueError, match='unambiguous ACGT'):
        project('AACCGG', 'NACCGG')


def test_all_retained_candidates_must_agree():
    def alias(target, p):
        return dict(target=target, sequence='AACCGG', mature_start=p)
    def hit(target, **kwargs):
        return dict(sequence_id='q', target=target, qstart=1, qend=6,
                    tstart=1, tend=6, qlen=6, tlen=6, qaln='AACCGG', taln='AACCGG', **kwargs)
    row = project(extra_contexts=[alias('alias', 4)], extra_hits=[hit('alias')])[0]
    assert row['status'] == 'resolved' and row['candidate_count'] == 2
    row = project(extra_contexts=[alias('other', 3)], extra_hits=[hit('other')])[0]
    assert row['status'] == 'conflicting_boundary' and row['mature_start'] is None
    unsupported = hit('t') | dict(qstart=4, tstart=4, qaln='CGG', taln='CGG')
    row = project(extra_hits=[unsupported])[0]
    assert row['status'] == 'unsupported_boundary'
    assert row['alignment_count'] == 2 and row['supported_alignment_count'] == 1


def test_invalid_sources_coordinates_and_schema_rejected():
    for kwargs in (dict(p=1), dict(p=True), dict(query='AA-CCGG'),
                   dict(qstart=2), dict(context='TTTTTT'), dict(mask=[]), dict(mask=[1]*6)):
        with pytest.raises(ValueError):
            project(**kwargs)
    with pytest.raises(ValueError, match='unknown query'):
        project(extra_hits=[dict(sequence_id='unknown', target='t')])
    with pytest.raises(ValueError, match='duplicate target'):
        project(extra_contexts=[dict(target='t', sequence='AACCGG', mature_start=4)])
    with pytest.raises(ValueError, match='coding-oriented'):
        project(extra_hits=[dict(sequence_id='q', target='t', qstart=6, qend=1, qlen=6)])
    assert project_transcript_boundaries([], [], []).shape == (0, 8)


def test_native_best_bits_retains_all_equal_alternatives_only():
    from vdjaqp.transcript_boundary import _project_best_boundaries
    queries = [dict(sequence_id='q', sequence='AACCGG', retained_primary_exact_match_mask=[True]*6)]
    contexts = [dict(target=t, sequence='AACCGG', mature_start=p) for t, p in [('a',4), ('b',4), ('c',3)]]
    high = dict(sequence_id='q', target='a', qstart=1, qend=6, tstart=1, tend=6,
                qlen=6, tlen=6, qaln='AACCGG', taln='AACCGG', bits=30.)
    weak = high | dict(target='b', qstart=4, tstart=4, qaln='CGG', taln='CGG', bits=10.)
    row = _project_best_boundaries(queries, contexts, [high, weak]).row(0, named=True)
    assert row['status']=='resolved' and row['mature_start']==4
    assert (row['returned_alignment_count'],row['alignment_count'],row['lower_score_alignment_count'])==(2,1,1)
    assert _project_best_boundaries(queries, contexts, [weak, high]).row(0,named=True)==row
    assert _project_best_boundaries(queries,contexts,[high,weak|dict(bits=30.)])['status'][0]=='unsupported_boundary'
    assert _project_best_boundaries(queries,contexts,[high,high|dict(target='c')])['status'][0]=='conflicting_boundary'
    assert _project_best_boundaries(queries,contexts,[high,high|dict(target='b')])['alignment_count'][0]==2
    assert _project_best_boundaries(queries,contexts,[])['status'][0]=='no_alignment'
