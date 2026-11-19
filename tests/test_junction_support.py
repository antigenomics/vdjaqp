"""Fixed contextual placements; no sequence search, annotation or Pgen call."""
import polars as pl
import pytest
from vdjaqp.junction_support import contextual_n_support,inferred_n_footprints,observed_junction_classes


def fixture():
    nt='TGT'+'GCC'*3+'CTG'+'TTT'+'TGG'
    c=pl.DataFrame([dict(clone_id='c',seq='QQCAAALFWEE',cdr3_start=2,cdr3_seq='CAAALFW')])
    a=pl.DataFrame([dict(clone_id='c',cdr3_nt=nt,v_end_nt=3,j_start_nt=18,
        d_start_nt=16,d_end_nt=18,np1=nt[3:15],np2=None)])
    scores=pl.DataFrame({'clone_id':['c'],'confirmed_nats':[20.]})
    support=pl.DataFrame(schema={'observation_id':pl.String,'clone_id':pl.String})
    return c,a,scores,support


def run(seq='CAAAIFW',minimum=3,**kwargs):
    c,a,scores,support=fixture()
    p=pl.DataFrame({'peptide_id':['p'],'seq':[seq]})
    placements=pl.DataFrame([dict(observation_id='o',peptide_id='p',clone_id='c',
        start=2,end=9,is_germline=False)])
    return contextual_n_support(c,p,placements,kwargs.get('annotations',a),accepted_scores=scores,
        accepted_support=kwargs.get('accepted_support',support),min_consecutive=minimum).row(0,named=True)


def test_contiguous_support_mismatch_break_and_il_equivalence():
    good=run();assert good['longest_n_run']==4 and good['rescue_eligible']
    broken=run('CAAVIFW');assert broken['agreeing_n_residues']==3 and broken['longest_n_run']==2
    assert not broken['rescue_eligible'] and run('CAAVIFW',2)['rescue_eligible']
    with pytest.raises(ValueError):run(minimum=0)


def test_mixed_d_and_missing_geometry_are_not_n_matches():
    c,a,_,_=fixture();r=a.row(0,named=True);r.update(v_end_nt=4,np1=r['cdr3_nt'][4:15])
    mixed=run(annotations=pl.DataFrame([r]));assert mixed['mixed_codon_residues']==1
    assert mixed['agreeing_n_residues']==3 and mixed['d_overlapping_residues']==1
    assert run(annotations=a.head(0))['geometry_status']=='missing_geometry'
    assert not run(annotations=a.head(0))['rescue_eligible']
    r['np1']='AAA'
    with pytest.raises(ValueError,match='slices'):inferred_n_footprints(c,pl.DataFrame([r]))


def test_rescue_mask_and_germline_competitors_remain_visible():
    c,a,scores,support=fixture()
    credited=pl.DataFrame({'observation_id':['o'],'clone_id':['c']})
    assert run(accepted_support=credited)['excluded_mapped_observation']
    assert not run(accepted_support=credited)['rescue_eligible']
    p=pl.DataFrame({'peptide_id':['p'],'seq':['CAAALFW']})
    rows=pl.DataFrame([dict(observation_id='o',peptide_id='p',clone_id='c',start=2,end=9,is_germline=v) for v in (False,True)])
    out=contextual_n_support(c,p,rows,a,accepted_scores=scores,accepted_support=support)
    assert out['rescue_eligible'].to_list()==[True,False] and out.height==2
    assert out['is_germline'].to_list()==[False,True]
    germline=c.with_columns(pl.lit(None).alias('cdr3_seq'),pl.lit(None).alias('cdr3_start'))
    unknown=contextual_n_support(germline,p,rows,a.head(0),accepted_scores=scores,accepted_support=support)
    assert unknown.height==2 and not unknown['rescue_eligible'].any()


def test_incomplete_or_discordant_observed_junctions_stay_unresolved():
    c,a,_,_=fixture()
    assert observed_junction_classes(c)=={'c':'CAAALFW'}
    for bad in [c.with_columns(pl.lit(0).alias('cdr3_start')),
                c.with_columns(pl.lit('AAALF').alias('cdr3_seq'))]:
        assert observed_junction_classes(bad)=={}
        assert inferred_n_footprints(bad,a.head(0))['c']['status']=='missing_geometry'
