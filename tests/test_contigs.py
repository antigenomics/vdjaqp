import pytest
from vdjaqp.contigs import adapt_contigs


def fixture(reverse=False):
    from arda.refbuild.translate import reverse_complement
    nt='GCT'*6+'TGTGCTTTT'+'GCT'*6
    mask=[True]*len(nt)
    submitted=reverse_complement(nt) if reverse else nt
    raw=[dict(sequence_id='s_c',sequence=submitted,scenario='s',contig_id='c')]
    ann=[dict(sequence_id='s_c',sequence=submitted,rev_comp='T' if reverse else 'F',locus='IGH',junction_aa='CAF',v_call='IGHV1-2*01',j_call='IGHJ4*01',fwr3_start=1,fwr3_end=21,fwr4_start=25,fwr4_end=45)]
    support={'s_c':dict(sequence=submitted,retained_primary_exact_match_mask=mask)}
    return raw,ann,support


def test_orientation_support_anchors_and_catalog():
    for reverse in (False,True):
        raw,ann,support=fixture(reverse)
        out=adapt_contigs(raw,ann,support)
        c=out['candidates'].row(0,named=True)
        assert c['seq']=='AAAAAACAF'+'AAAAAA'
        assert c['region_ids'][6:9]==[4,5,7]
        assert c['cdr3_start']==6 and c['cdr3_end']==9
        assert out['catalog'].height==1 and out['search_entries'].height==1
        assert not {'family_id','chain_id','antibody_id'} & set(c)
    raw,ann,support=fixture()
    support['s_c']['retained_primary_exact_match_mask'][19]=False
    out=adapt_contigs(raw,ann,support)
    assert out['candidates'].is_empty()
    assert out['catalog']['status'][0]=='no_supported_complete_junction'
    ann[0]['locus']='IGL'
    out=adapt_contigs(raw,ann,support)
    assert out['catalog'].height==1 and out['candidates'].is_empty()
    assert out['search_entries'].height==2
    support['s_c']['retained_primary_exact_match_mask']=[]
    with pytest.raises(ValueError):adapt_contigs(raw,ann,support)

def test_asymmetric_reverse_mask_and_partial_terminal_codon():
    for reverse in (False,True):
        raw,ann,support=fixture(reverse)
        index=len(raw[0]['sequence'])-20 if reverse else 19
        support['s_c']['retained_primary_exact_match_mask'][index]=False
        assert adapt_contigs(raw,ann,support)['candidates'].is_empty()
    raw,ann,support=fixture()
    for r in (raw[0],ann[0],support['s_c']):r['sequence']+='A'
    support['s_c']['retained_primary_exact_match_mask'].append(False)
    out=adapt_contigs(raw,ann,support)
    assert len(out['candidates']['seq'][0])==15
    assert out['catalog']['raw_nt_length'][0]==46


@pytest.mark.parametrize('orientation', ['', None])
def test_unannotated_contig_retained_without_coding_orientation(orientation):
    # Native unannotated AIRR rows contain the submitted sequence and empty fields.
    nt='ACGTACGTACGT'
    raw=[dict(sequence_id='unannotated',sequence=nt)]
    ann=[dict(sequence_id='unannotated',sequence=nt,rev_comp=orientation,
              locus='',junction_aa='',v_call='',j_call='')]
    support={'unannotated':dict(sequence=nt,retained_primary_exact_match_mask=[False]*len(nt))}
    out=adapt_contigs(raw,ann,support)
    row=out['catalog'].row(0,named=True)
    assert row['status']=='unavailable_coding_orientation'
    assert row['raw_nt_length']==len(nt) and row['locus']=='unknown'
    assert row['rev_comp'] is row['coding_frame'] is row['coding_aa_length'] is None
    assert row['supported_annotated_aa']==0
    assert out['candidates'].is_empty() and out['search_entries'].is_empty()
    assert raw[0]['sequence']==ann[0]['sequence']==support['unannotated']['sequence']==nt
    support['unannotated']['retained_primary_exact_match_mask']=[]
    with pytest.raises(ValueError,match='invalid retained-primary support mask'):
        adapt_contigs(raw,ann,support)
    support['unannotated']['retained_primary_exact_match_mask']=[False]*len(nt)
    ann[0]['sequence']='ACGT'
    with pytest.raises(ValueError,match='sequence source mismatch'):
        adapt_contigs(raw,ann,support)


@pytest.mark.parametrize('orientation', ['unknown', 'N', ' '])
def test_invalid_nonblank_coding_orientation_rejected(orientation):
    raw,ann,support=fixture()
    ann[0]['rev_comp']=orientation
    with pytest.raises(ValueError,match='missing coding orientation'):
        adapt_contigs(raw,ann,support)
