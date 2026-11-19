"""J–C spanning peptides give protein-level isotype evidence on toy candidates."""
import polars as pl
import pytest

from vdjaqp.isotype_evidence import (IGHC_CH1, c_name_isotype, isotype_evidence,
                                    junction_constant_placements)

FR4 = 'WGQGTLVTVSS'


def candidate(clone_id, cdr3, constant, c_name):
    seq = 'EVQL' + cdr3 + FR4[1:] + constant
    labels = [4] * 4 + [4] + [5] * (len(cdr3) - 2) + [7] + [7] * (len(FR4) - 1) + [8] * len(constant)
    return dict(clone_id=clone_id, seq=seq, cdr3_seq=cdr3, region_ids=labels, c_name=c_name)


def toy():
    cands = pl.DataFrame([
        candidate('g1', 'CARDYW', 'ASTKGPSVFPLAPSSK', 'IGHG1'),
        candidate('g2', 'CTRIIW', 'ASTKGPSVFPLAPCSR', 'IGHG2'),
        candidate('a1', 'CAKEEW', 'ASPTSPKVFPLSLCST', 'IGHA1'),
        candidate('m_bad', 'CVRGGW', 'GSASAPTLFPLV', 'IGHG1'),
        candidate('v_only', 'CASSSW', '', 'unknown'),
    ])
    return cands


def place(cands, clone_id, peptide, peptide_id=None):
    seq = cands.filter(pl.col('clone_id') == clone_id)['seq'][0]
    start = seq.index(peptide)
    return dict(clone_id=clone_id, peptide_id=peptide_id or peptide, start=start, end=start + len(peptide))


def test_spanning_requires_fr4_and_native_constant():
    cands = toy()
    rows = [place(cands, 'g1', 'TVSSASTKGPSVFPLAPSSK'),  # spans, subclass-specific S at CH1 14
            place(cands, 'g1', 'GTLVTVSS'),               # FR4 only
            place(cands, 'g1', 'STKGPSVF'),                # constant only
            place(cands, 'v_only', 'GTLVTVSS')]
    jc = junction_constant_placements(pl.DataFrame(rows), cands)
    assert jc.height == 1
    row = jc.row(0, named=True)
    assert row['c_offset'] == 0 and row['c_residues'] == 'ASTKGPSVFPLAPSSK'
    assert row['consistent_genes'] == 'IGHG1' and row['consistent_isotypes'] == 'IgG'


def test_shared_igg_prefix_is_isotype_level_and_subclass_when_resolved():
    cands = toy()
    rows = [place(cands, 'g1', 'VTVSSASTKGPS', 'p_shared'),
            place(cands, 'g2', 'VTVSSASTKGPS', 'p_shared'),
            place(cands, 'g2', 'SSASTKGPSVFPLAPCSR', 'p_g2'),
            place(cands, 'a1', 'TVSSASPTSPK', 'p_a')]
    out = isotype_evidence(pl.DataFrame(rows), cands)
    by = {r['identity_class']: r for r in out.iter_rows(named=True)}
    g1, g2, a1 = by['CARDYW'], by['CTRLLW'], by['CAKEEW']
    assert g1['protein_isotype'] == 'IgG' and g1['protein_subclass'] is None
    assert g1['jc_peptides'] == 1 and g1['jc_peptides_unique_to_class'] == 0
    assert g1['identity_linked_isotype'] is None
    # CH1 residues 1-16 ending in CSR are shared by IGHG2, IGHG3 and IGHG4.
    assert g2['protein_isotype'] == 'IgG' and g2['protein_subclass'] is None
    assert g2['jc_peptides'] == 2 and g2['jc_peptides_unique_to_class'] == 1
    assert g2['identity_linked_isotype'] == 'IgG'
    assert a1['protein_isotype'] == 'IgA' and a1['c_name_agreement'] == 'agree'
    # A J–C placement on one class is not class-unique if the peptide also sits elsewhere.
    elsewhere = rows + [dict(clone_id='v_only', peptide_id='p_a', start=0, end=4)]
    a1b = isotype_evidence(pl.DataFrame(elsewhere), cands).filter(pl.col('identity_class') == 'CAKEEW').row(0, named=True)
    assert a1b['jc_peptides_unique_to_class'] == 0 and a1b['identity_linked_isotype'] is None
    assert a1b['protein_isotype'] == 'IgA'
    assert by['CASSSW']['jc_peptides'] == 0 and by['CASSSW']['protein_isotype'] is None
    assert by['CASSSW']['c_name_agreement'] == 'not_called'
    assert out['accepted'].null_count() == out.height


def test_disagreement_with_c_name_and_reference_mismatch():
    cands = toy()
    rows = [place(cands, 'm_bad', 'TVSSGSASAPTL'),
            dict(clone_id='g1', peptide_id='mut', start=place(cands, 'g1', 'VTVSSASTK')['start'],
                 end=place(cands, 'g1', 'VTVSSASTK')['end'], peptide='VTVSSAWTK')]
    out = isotype_evidence(pl.DataFrame(rows), cands)
    by = {r['identity_class']: r for r in out.iter_rows(named=True)}
    assert by['CVRGGW']['protein_isotype'] == 'IgM'
    assert by['CVRGGW']['c_name_isotype'] == 'IgG' and by['CVRGGW']['c_name_agreement'] == 'disagree'
    assert by['CARDYW']['jc_peptides'] == 1 and by['CARDYW']['reference_consistent_peptides'] == 0
    assert by['CARDYW']['protein_isotype'] is None


def test_unresolved_conflicting_and_identity_rule():
    cands = toy()
    # One constant residue 'A' fits IgG, IgA, IgD and IgE alike.
    rows = [place(cands, 'g1', 'GTLVTVSSA'),
            place(cands, 'a1', 'TVSSASPTSPK')]
    extra = candidate('a_twin', 'CARDYW', 'ASPTSPKVFPL', 'IGHA1')
    cands2 = pl.concat([cands, pl.DataFrame([extra], schema=cands.schema)])
    rows2 = rows + [place(cands2, 'a_twin', 'VTVSSASPTS'), place(cands2, 'g1', 'VTVSSASTKG')]
    scores = pl.DataFrame({'clone_id': ['g1', 'a1', 'a_twin'], 'confirmed_nats': [25.0, 3.0, 1.0]})
    out = isotype_evidence(pl.DataFrame(rows), cands, scores)
    by = {r['identity_class']: r for r in out.iter_rows(named=True)}
    assert by['CARDYW']['protein_isotype'] == 'unresolved' and by['CARDYW']['accepted'] is True
    assert by['CAKEEW']['accepted'] is False and by['CASSSW']['accepted'] is False
    out2 = isotype_evidence(pl.DataFrame(rows2), cands2, scores, threshold=0.5)
    twin = out2.filter(pl.col('identity_class') == 'CARDYW').row(0, named=True)
    assert twin['records'] == 2 and twin['protein_isotype'] == 'conflicting'
    assert twin['c_name_agreement'] == 'not_called' and twin['c_name_isotype'] is None


def test_identity_class_column_and_validation():
    cands = toy().with_columns(pl.col('clone_id').alias('identity_class'))
    out = isotype_evidence(pl.DataFrame([place(cands, 'g1', 'TVSSASTKGP')]), cands)
    assert set(out['identity_class']) == set(cands['clone_id'])
    bad = pl.DataFrame([dict(clone_id='g1', peptide_id='x', start=30, end=40)])
    with pytest.raises(ValueError):
        junction_constant_placements(bad, toy())


def test_c_name_parsing_and_reference_prefixes():
    assert c_name_isotype('IGHG1*01') == 'IgG'
    assert c_name_isotype('IGHG1,IGHG3') == 'IgG'
    assert c_name_isotype('IGHG1,IGHM') is None
    assert c_name_isotype('unknown') is None and c_name_isotype(None) is None
    assert {IGHC_CH1[g][:13] for g in ('IGHG1', 'IGHG2', 'IGHG3', 'IGHG4')} == {'ASTKGPSVFPLAP'}
    assert IGHC_CH1['IGHM'].startswith('GSASAPTLFP') and IGHC_CH1['IGHA1'].startswith('ASPTSPKVFPL')
