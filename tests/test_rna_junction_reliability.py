"""Optional RNA-class calibration and conservative penalty; artificial data only."""
import numpy as np
import polars as pl
import pytest


def candidates():
    frame=pl.DataFrame({'clone_id':['a','a_alias','b','c'],
        'parent_clone_id':['p1','p1','p2','p3'],'cdr3_seq':['CAAW','CAAW','CCDW','CEFW'],
        'representative_policy':['native_root']*4,'native_root_matching_observations':[1,1,8,0]})
    identity=frame.select('clone_id',pl.col('cdr3_seq').alias('identity_class'))
    return frame,identity


def scores():
    return pl.DataFrame({'identity_class':['CAAW','CCDW','CEFW'],
        'profiled_removal_lower_bound':[2.,2.,None],
        'profile_score_available':[True,True,False],
        'identifiable':[True,True,False],
        'profile_status':['profiled_certified_bound','profiled_certified_bound','observationally_equivalent']})


def test_ties_smoothed_once_per_bin_with_class_weights_and_permutation_invariance():
    from vdjaqp.errormodel import fit_rna_junction_reliability
    x=np.array([1,1,8]);g=np.array([True,False,True])
    knots,q=fit_rna_junction_reliability(x,g)
    assert np.array_equal(knots,[1.,8.]) and np.allclose(q,[.5,.75])
    for order in [[2,1,0],[1,0,2]]:
        other,p=fit_rna_junction_reliability(x[order],g[order])
        assert np.array_equal(knots,other) and np.array_equal(q,p)
    # Inverted bins pool by their class counts (3 versus 1), not equal bin weights.
    knots,q=fit_rna_junction_reliability([1,1,1,8],[True,True,True,False])
    expected=1.-(3*(.5/4)+1*(1.5/2))/4
    assert np.allclose(q,expected)


@pytest.mark.parametrize('labels',[[True]*4,[False]*4])
def test_smoothing_keeps_reliability_strictly_between_zero_and_one(labels):
    from vdjaqp.errormodel import fit_rna_junction_reliability
    x,q=fit_rna_junction_reliability([0,1,8,100],labels)
    assert np.isfinite(q).all() and ((q>0)&(q<1)).all()
    assert (np.diff(q)>=0).all()


@pytest.mark.parametrize('support,genuine',[
    ([],[]),([-1],[True]),([None],[True]),([np.nan],[True]),([1.5],[True]),
    ([True],[True]),([1],[None]),([1],[2]),([1,2],[True]),([[1]],[True])])
def test_invalid_calibration_inputs_fail_closed(support,genuine):
    from vdjaqp.errormodel import fit_rna_junction_reliability
    with pytest.raises(ValueError):fit_rna_junction_reliability(support,genuine)


def test_full_class_support_maximum_is_alias_and_row_invariant():
    from vdjaqp.aggregate import identity_rna_support
    frame,identity=candidates()
    result=identity_rna_support(frame,identity)
    other=identity_rna_support(frame.filter(pl.col('clone_id')!='a_alias'),identity.filter(pl.col('clone_id')!='a_alias'))
    assert result.equals(other) and result.equals(identity_rna_support(frame.reverse(),identity.reverse()))
    assert result['rna_root_observations'].to_list()==[1,8,0]
    assert result['native_parents'].to_list()==[1,1,1]
    more=frame.filter(pl.col('clone_id')=='a').with_columns(pl.lit('another').alias('clone_id'),
        pl.lit('p4').alias('parent_clone_id'),pl.lit(3,dtype=pl.Int64).alias('native_root_matching_observations'))
    more_ids=more.select('clone_id',pl.col('cdr3_seq').alias('identity_class'))
    combined=identity_rna_support(pl.concat([frame,more]),pl.concat([identity,more_ids]))
    assert combined['rna_root_observations'][0]==3 and combined['native_parents'][0]==2


def test_penalty_is_bounded_truth_free_and_retains_unavailable_ambiguity_population():
    from vdjaqp.aggregate import identity_rna_support,apply_rna_junction_reliability
    frame,identity=candidates();support=identity_rna_support(frame,identity)
    before=scores()
    result=apply_rna_junction_reliability(before,support,[0.,1.,8.],[.25,.5,.9])
    assert result.select(before.columns).equals(before)
    assert result.height==3
    assert np.allclose(result['rna_junction_reliability'],[.5,.9,.25])
    assert np.allclose(result['rna_adjusted_removal_lower_bound'][:2],[2+np.log(.5),2+np.log(.9)])
    assert result['rna_adjusted_removal_lower_bound'][2] is None
    assert result['rna_adjusted_score_available'].to_list()==[True,True,False]
    assert (result['rna_adjusted_removal_lower_bound'][:2]<=before['profiled_removal_lower_bound'][:2]).all()
    endpoints=apply_rna_junction_reliability(before,support,[2.,5.],[.4,.8])
    assert endpoints['rna_junction_reliability'].to_list()==[.4,.8,.4]
    parity=apply_rna_junction_reliability(before,support,[0.,8.],[1.,1.])
    assert parity['rna_adjusted_removal_lower_bound'].equals(before['profiled_removal_lower_bound'])


def test_artificial_false_rna_class_penalized_without_dropping_true_class():
    from vdjaqp.errormodel import fit_rna_junction_reliability
    from vdjaqp.aggregate import identity_rna_support,apply_rna_junction_reliability
    knots,q=fit_rna_junction_reliability([1]*3+[8]*3,[False]*3+[True]*3)
    frame,identity=candidates()
    result=apply_rna_junction_reliability(scores(),identity_rna_support(frame,identity),knots,q)
    assert result['rna_adjusted_removal_lower_bound'][0]==0.
    assert result['rna_adjusted_removal_lower_bound'][1]>1.
    assert result.height==3


@pytest.mark.parametrize('change',['negative','missing','float','policy','parent','mapping','alias_conflict'])
def test_bad_native_support_metadata_fails_closed(change):
    from vdjaqp.aggregate import identity_rna_support
    frame,identity=candidates()
    if change=='negative':frame=frame.with_columns(pl.lit(-1).alias('native_root_matching_observations'))
    if change=='missing':frame=frame.with_columns(pl.lit(None,dtype=pl.Int64).alias('native_root_matching_observations'))
    if change=='float':frame=frame.with_columns(pl.col('native_root_matching_observations').cast(pl.Float64))
    if change=='policy':frame=frame.with_columns(pl.lit('longest').alias('representative_policy'))
    if change=='parent':frame=frame.with_columns(pl.lit(None,dtype=pl.String).alias('parent_clone_id'))
    if change=='mapping':identity=identity.filter(pl.col('clone_id')!='c')
    if change=='alias_conflict':frame=frame.with_columns(pl.when(pl.col('clone_id')=='a_alias').then(3).otherwise(pl.col('native_root_matching_observations')).alias('native_root_matching_observations'))
    with pytest.raises(ValueError):identity_rna_support(frame,identity)


@pytest.mark.parametrize('knots,q',[
    ([0,8],[0.,.9]),([0,8],[.5,1.1]),([0,8],[.9,.5]),([8,0],[.5,.9]),
    ([0,0],[.5,.9]),([0,8],[.5,np.nan]),([],[]),([0],[.5,.9])])
def test_invalid_frozen_curve_fails_closed(knots,q):
    from vdjaqp.aggregate import identity_rna_support,apply_rna_junction_reliability
    frame,identity=candidates()
    with pytest.raises(ValueError):apply_rna_junction_reliability(scores(),identity_rna_support(frame,identity),knots,q)


@pytest.mark.parametrize('change',['missing_class','wrong_units','missing_count','negative_count'])
def test_application_requires_whole_population_and_valid_support(change):
    from vdjaqp.aggregate import identity_rna_support,apply_rna_junction_reliability
    frame,identity=candidates();support=identity_rna_support(frame,identity)
    if change=='missing_class':support=support.head(2)
    if change=='wrong_units':support=support.with_columns(pl.lit('molecules').alias('rna_support_units'))
    if change=='missing_count':support=support.with_columns(pl.lit(None,dtype=pl.Int64).alias('rna_root_observations'))
    if change=='negative_count':support=support.with_columns(pl.lit(-1).alias('rna_root_observations'))
    with pytest.raises(ValueError):apply_rna_junction_reliability(scores(),support,[0,8],[.5,.9])


def test_application_preserves_score_order_independent_of_support_order():
    from vdjaqp.aggregate import identity_rna_support,apply_rna_junction_reliability
    frame,identity=candidates();support=identity_rna_support(frame,identity)
    first=apply_rna_junction_reliability(scores(),support,[0,8],[.5,.9])
    second=apply_rna_junction_reliability(scores().reverse(),support.reverse(),[0,8],[.5,.9])
    assert first.equals(second.reverse())
