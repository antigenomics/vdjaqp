import pytest
from vdjaqp.transcript_context import select_contexts, leader_rna_template


def test_donor_hierarchy_ties_no_hit_and_unchanged_mature_payload():
    sources = [dict(source_id='s1', locus='IGH', core_nt='AAA', status='valid',
                    prefix_nt='ATGGCT', leader_nt='ATGGCT'),
               dict(source_id='s2', locus='IGH', core_nt='CCC', status='valid',
                    prefix_nt='ATGGCC', leader_nt='ATGGCC')]
    def native(name, gene, core, ids):
        return dict(native_allele=name, gene=gene, core_nt=core, locus='IGH', exact_sources=ids)
    natives = [native('exact', 'g1', 'AAA', ['s1']), native('exact2', 'g2', 'CCC', ['s2']),
               native('group', 'g1', 'AAT', []), native('nearest', 'g3', 'ACC', []),
               native('missing', 'g4', 'TTT', [])]
    hits = [dict(native_allele='nearest', source_id=s, bits=12.0) for s in ('s1', 's2')]
    rows, provenance = select_contexts(natives, sources, hits)
    assert [r['method'] for r in rows] == ['STRICT_EXACT', 'STRICT_EXACT', 'IMPUTED_GROUP', 'IMPUTED_ALIGNMENT', 'UNSUPPORTED']
    assert rows[3]['confidence'] == 'ambiguous_tied_prefixes'
    assert len([r for r in provenance if r['native_allele'] == 'nearest']) == 2
    assert rows == select_contexts(natives, list(reversed(sources)), list(reversed(hits)))[0]
    assert rows[4]['prefix_nt'] is None
    result = leader_rna_template('CAG', 'GCTTAA', 'QA', rows[0], v_trim_5=0)
    assert result['rnaseq_nt'] == 'ATGGCTCAGGCTTAA' and result['mature_nt_start'] == 7
    for variable, constant, protein, trim in [('CAG', 'GCTTAA', 'QA', 1), ('CAA', 'GCTTAA', 'wrong', 0)]:
        with pytest.raises(ValueError):
            leader_rna_template(variable, constant, protein, rows[0], v_trim_5=trim)
    with pytest.raises(ValueError, match='cross-locus'):
        select_contexts([dict(natives[3], locus='IGK'), *natives[:2]], sources, hits)
    for malformed in (dict(natives[0], exact_sources=['unknown']),
                      dict(natives[0], core_nt=''), dict(natives[0], locus='unknown'),
                      dict(natives[0], native_allele=''), dict(natives[0], gene=None)):
        with pytest.raises(ValueError): select_contexts([malformed], sources)
    with pytest.raises(ValueError): select_contexts(natives, [dict(sources[0], status=None)])
    with pytest.raises(ValueError): select_contexts(natives, [dict(sources[0], status='unknown')])


def test_leader_rna_truth_validation_preserves_mature_regions():
    from test_simulation_truth import fixture
    from vdjaqp.simulation_truth import regions_from_annotations
    chain, ann = fixture()
    original = regions_from_annotations([chain], [ann])
    mature = chain['variable_nt']+'CCTCCTTAA'
    chain.update(rna_model='leader_imputed_v1', rna_prefix_nt='ATGGCT', rna_leader_nt='ATGGCT',
                 rna_mature_nt_start=7, rna_leader_nt_start=1, rnaseq_nt='ATGGCT'+mature)
    result = regions_from_annotations([chain], [ann])
    assert original[0].equals(result[0]) and original[1].equals(result[1])
    chain['rna_mature_nt_start'] = 4
    with pytest.raises(ValueError, match='leader RNA'):
        regions_from_annotations([chain], [ann])


@pytest.mark.parametrize('reverse', [False, True])
def test_observed_boundary_preserves_unaligned_read_supported_fr1(reverse):
    import polars as pl
    from test_contigs import fixture
    from vdjaqp.contigs import adapt_contigs
    from arda.refbuild.translate import reverse_complement
    raw, ann, supports = fixture()
    name = raw[0]['sequence_id']
    nt = 'ATG'*3 + raw[0]['sequence']
    raw[0]['sequence'] = reverse_complement(nt) if reverse else nt
    ann[0].update(sequence=raw[0]['sequence'], rev_comp='T' if reverse else 'F',
                  fwr1_start=16, fwr1_end=18, fwr3_start=16, fwr3_end=30,
                  fwr4_start=34, fwr4_end=54)
    supports[name].update(sequence=raw[0]['sequence'], retained_primary_exact_match_mask=[True]*len(nt))
    ordinary = adapt_contigs(raw, ann, supports)['candidates'].row(0, named=True)
    boundaries = pl.DataFrame([dict(sequence_id=name, mature_start=10, status='resolved')])
    candidate = adapt_contigs(raw, ann, supports, boundaries=boundaries)['candidates'].row(0, named=True)
    assert candidate['seq'] == 'AAAAAACAF'+'AAAAAA'
    assert candidate['coding_aa_start'] == 3 and ordinary['coding_aa_start'] == 5
    assert candidate['source_nt_start'] == (0 if reverse else 9)
    boundaries = boundaries.with_columns(pl.lit('unsupported_boundary').alias('status'))
    assert adapt_contigs(raw, ann, supports, boundaries=boundaries)['candidates'].row(0, named=True) == ordinary
