"""Experimental source-class competition retains aliases, null and ambiguity."""
import polars as pl
import numpy as np
from vdjaqp.aggregate import compete_identity_classes


def test_aliases_cannot_multiply_source_evidence_and_ties_unresolved():
    identities=pl.DataFrame({'clone_id':['a','alias','b'],'identity_class':['A','A','B']})
    peptides=pl.DataFrame({'peptide_id':['p'],'spectrum_id':['s'],'alc':[1.], 'alc_rank':[1.], 'n_reads':[1]})
    scored=pl.DataFrame({'clone_id':['a','alias','b'],'peptide_id':['p']*3,'llr':[10.,10.,10.]})
    out=compete_identity_classes(scored,peptides,identities)
    assert out.height==2 and not out['unique_map_above_null'].any()
    assert out['responsibility'][0]==out['responsibility'][1]
    assert out.equals(compete_identity_classes(scored.filter(pl.col('clone_id')!='alias'),peptides,identities))


def test_null_rejects_unique_weak_class():
    identities=pl.DataFrame({'clone_id':['a'],'identity_class':['A']})
    peptides=pl.DataFrame({'peptide_id':['p'],'spectrum_id':['s'],'alc':[1.], 'alc_rank':[1.], 'n_reads':[1]})
    scored=pl.DataFrame({'clone_id':['a'],'peptide_id':['p'],'llr':[-1.]})
    out=compete_identity_classes(scored,peptides,identities)
    assert not out['unique_map_above_null'].any()
    assert out['responsibility'][0]<.5


def test_spectrum_confidence_and_source_competition_are_separate():
    identities=pl.DataFrame({'clone_id':['a','b'],'identity_class':['A','B']})
    peptides=pl.DataFrame({'peptide_id':['p','q'],'spectrum_id':['s','t'],
        'alc':[1.,.1], 'alc_rank':[1.,1.], 'n_reads':[1,1]})
    scored=pl.DataFrame({'clone_id':['a','b','a','b'],'peptide_id':['p','p','q','q'],
        'llr':[2.,1.,2.,1.]})
    out=compete_identity_classes(scored,peptides,identities)
    assert out.filter(pl.col('spectrum_id')=='s')['unique_map_above_null'].sum()==1
    assert not out.filter(pl.col('spectrum_id')=='t')['unique_map_above_null'].any()
    assert out.group_by('spectrum_id').agg(pl.col('responsibility').sum())['responsibility'].max()<1


def test_incomplete_identity_mapping_fails_closed():
    import pytest
    identities=pl.DataFrame({'clone_id':['a'],'identity_class':['A']})
    peptides=pl.DataFrame({'peptide_id':['p'],'spectrum_id':['s'],'alc':[1.], 'alc_rank':[1.], 'n_reads':[1]})
    scored=pl.DataFrame({'clone_id':['missing'],'peptide_id':['p'],'llr':[1.]})
    with pytest.raises(ValueError,match='lacks identity'):
        compete_identity_classes(scored,peptides,identities)


def test_mixture_retains_two_sources_and_rejects_shared_decoy():
    from vdjaqp.aggregate import fit_identity_mixture
    peptides=pl.DataFrame({'peptide_id':['p','q','r'],'spectrum_id':['s','t','u'],
        'seq':['ACDEFGHIK','YYYYYYYYY','RRRRRRRRR']})
    competition=pl.DataFrame({'spectrum_id':['s','s','t','t','u','u','u'],
        'identity_class':['A','D','B','D','A','B','D'],
        'log_weight':[10.,0.,8.,0.,3.,3.,0.]})
    out,receipt=fit_identity_mixture(competition,peptides)
    assert receipt['converged']
    assert (out.filter(pl.col('identity_class').is_in(['A','B']))['mixture_weight']>.1).all()
    assert out.filter(pl.col('identity_class')=='D')['mixture_weight'][0]<1e-6
    assert np.isfinite(out['fixed_other_removal_contrast'].to_numpy()).all()
    other,_=fit_identity_mixture(competition.reverse(),peptides.reverse())
    assert np.allclose(out['mixture_weight'],other['mixture_weight'])


def test_single_dominant_source_preserves_log_null_and_finite_removal():
    from vdjaqp.aggregate import fit_identity_mixture
    peptides=pl.DataFrame({'peptide_id':['p'],'spectrum_id':['s'],'seq':['ACDEFGHIK']})
    competition=pl.DataFrame({'spectrum_id':['s'],'identity_class':['A'],'log_weight':[1000.]})
    out,receipt=fit_identity_mixture(competition,peptides)
    assert receipt['converged']
    assert out['mixture_weight'][0]==1.
    assert out['removal_contrast_available'][0]
    assert np.isclose(out['fixed_other_removal_contrast'][0],1000.)
    assert receipt['unavailable_contrasts']==0
    assert np.isfinite(receipt['log_null_weight'])


def test_huge_two_source_likelihood_has_stable_removal_and_monotone_fit():
    from vdjaqp.aggregate import fit_identity_mixture
    peptides=pl.DataFrame({'peptide_id':['p','q'],'spectrum_id':['s','t'],
        'seq':['ACDEFGHIK','YYYYYYYYY']})
    competition=pl.DataFrame({'spectrum_id':['s','t'],'identity_class':['A','B'],
        'log_weight':[1000.,1000.]})
    out,receipt=fit_identity_mixture(competition,peptides)
    assert receipt['converged']
    assert np.allclose(out['mixture_weight'],.5)
    assert out['removal_contrast_available'].all()
    assert np.allclose(out['fixed_other_removal_contrast'],1000.-3*np.log(2)-receipt['log_null_weight'])
    assert (np.diff(receipt['objective'])>=-1e-12).all()


def test_negative_likelihoods_fit_null_and_empty_classes_are_supported():
    from vdjaqp.aggregate import fit_identity_mixture
    peptides=pl.DataFrame({'peptide_id':['p'],'spectrum_id':['s'],'seq':['ACDEFGHIK']})
    competition=pl.DataFrame({'spectrum_id':['s'],'identity_class':['A'],'log_weight':[-10.]})
    out,receipt=fit_identity_mixture(competition,peptides)
    assert receipt['converged'] and receipt['null_weight']>1.-1e-8
    assert out['mixture_weight'][0]<1e-8
    empty=competition.clear()
    out,receipt=fit_identity_mixture(empty,peptides)
    assert out.height==0 and receipt['null_weight']==1.


def test_reacquired_il_equivalent_peptide_does_not_multiply_presence_evidence():
    from vdjaqp.aggregate import fit_identity_mixture
    peptides=pl.DataFrame({'peptide_id':['p','q'],'spectrum_id':['s','t'],
        'seq':['ACDEFGHIK','YYYYYYYYY']})
    competition=pl.DataFrame({'spectrum_id':['s','t'],'identity_class':['A','B'],
        'log_weight':[5.,2.]})
    baseline,_=fit_identity_mixture(competition,peptides)
    repeated=pl.concat([peptides,pl.DataFrame({'peptide_id':['r'],'spectrum_id':['u'],
        'seq':['ACDEFGHLK']})])
    edges=pl.concat([competition,pl.DataFrame({'spectrum_id':['u'],
        'identity_class':['A'],'log_weight':[5.]})])
    out,receipt=fit_identity_mixture(edges,repeated)
    assert receipt['weighted_observations']==2.
    assert np.allclose(out['mixture_weight'],baseline['mixture_weight'])
    assert np.allclose(out['fixed_other_removal_contrast'],baseline['fixed_other_removal_contrast'])


def test_weak_private_peptide_retains_second_true_component():
    from vdjaqp.aggregate import fit_identity_mixture
    peptides=pl.DataFrame({'peptide_id':['p','q'],'spectrum_id':['s','t'],
        'seq':['ACDEFGHIK','YYYYYYYYY']})
    competition=pl.DataFrame({'spectrum_id':['s','t','s','t'],
        'identity_class':['A','B','D','D'],'log_weight':[10.,1.,0.,0.]})
    out,receipt=fit_identity_mixture(competition,peptides)
    assert receipt['converged']
    true=out.filter(pl.col('identity_class').is_in(['A','B']))
    assert (true['mixture_weight']>.1).all()
    assert (true['fixed_other_removal_contrast']>0.).all()
    assert out.filter(pl.col('identity_class')=='D')['mixture_weight'][0]<1e-6


def test_sparse_removal_matches_explicit_all_observation_likelihood():
    from vdjaqp.aggregate import fit_identity_mixture
    peptides=pl.DataFrame({'peptide_id':['p','q','r'],'spectrum_id':['s','t','u'],
        'seq':['ACDEFGHIK','YYYYYYYYY','RRRRRRRRR']})
    competition=pl.DataFrame({'spectrum_id':['s','s','t','u'],
        'identity_class':['A','B','B','A'],'log_weight':[3.,1.,4.,-1.]})
    out,receipt=fit_identity_mixture(competition,peptides)
    weights=out['mixture_weight'].to_numpy()
    likelihood=np.array([[np.exp(3.),np.exp(1.)],[0.,np.exp(4.)],[np.exp(-1.),0.]])
    full=likelihood@weights+receipt['null_weight']
    expected=[]
    for i in range(2):
        keep=np.arange(2)!=i
        minus=(likelihood[:,keep]@weights[keep]+receipt['null_weight'])/(1.-weights[i])
        expected.append(np.log(full/minus).sum())
    assert np.allclose(out['fixed_other_removal_contrast'],expected)


def test_null_class_likelihood_fails_closed():
    import pytest
    identities=pl.DataFrame({'clone_id':['a'],'identity_class':['A']})
    peptides=pl.DataFrame({'peptide_id':['p'],'spectrum_id':['s'],'alc':[1.],
        'alc_rank':[1.],'n_reads':[1]})
    scored=pl.DataFrame({'clone_id':['a'],'peptide_id':['p'],
        'llr':pl.Series([None],dtype=pl.Float64)})
    with pytest.raises(ValueError,match='nonfinite'):
        compete_identity_classes(scored,peptides,identities)


def test_equal_likelihood_profiles_are_explicitly_unidentifiable():
    from vdjaqp.aggregate import fit_identity_mixture
    peptides=pl.DataFrame({'peptide_id':['p','q'],'spectrum_id':['s','t'],
        'seq':['ACDEFGHIK','YYYYYYYYY']})
    competition=pl.DataFrame({'spectrum_id':['s','t','s','t'],
        'identity_class':['A','A','B','B'],'log_weight':[5.,2.,5.,2.]})
    out,receipt=fit_identity_mixture(competition,peptides)
    assert receipt['converged']
    assert not out['identifiable'].any()
    assert out['identity_equivalence_group'].n_unique()==1
    assert out['identity_equivalence_size'].to_list()==[2,2]
    assert receipt['individually_unidentifiable_classes']==2
    shuffled,_=fit_identity_mixture(competition.reverse(),peptides.reverse())
    assert out.equals(shuffled)


def test_one_weak_distinguishing_observation_breaks_profile_equivalence():
    from vdjaqp.aggregate import fit_identity_mixture
    peptides=pl.DataFrame({'peptide_id':['p','q'],'spectrum_id':['s','t'],
        'seq':['ACDEFGHIK','YYYYYYYYY']})
    competition=pl.DataFrame({'spectrum_id':['s','s','t'],
        'identity_class':['A','B','A'],'log_weight':[5.,5.,.1]})
    out,receipt=fit_identity_mixture(competition,peptides)
    assert receipt['converged']
    assert out['identifiable'].all()
    assert out['identity_equivalence_group'].n_unique()==2
    assert out['identity_equivalence_size'].to_list()==[1,1]


def test_distinct_profiles_can_still_have_nonidentifiable_joint_explanations():
    """Separate V/J witnesses cannot phase the two alternative source pairings.

    This guards the documented limitation: exact-profile identifiability is only
    necessary, and fixed-other removal is not refitted class-presence evidence.
    """
    from vdjaqp.aggregate import fit_identity_mixture
    peptides=pl.DataFrame({'peptide_id':['p','q','r','t'],
        'spectrum_id':['V1','V2','J1','J2'],
        'seq':['ACDEFGHIK','YYYYYYYYY','RRRRRRRRR','PPPPPPPPP']})
    competition=pl.DataFrame({'spectrum_id':['V1','J1','V2','J2','V1','J2','V2','J1'],
        'identity_class':['A','A','B','B','C','C','D','D'],
        'log_weight':[5.]*8})
    out,receipt=fit_identity_mixture(competition,peptides)
    assert receipt['converged'] and out['identifiable'].all()
    assert (out['fixed_other_removal_contrast']>0.2).all()
    # A+B and C+D give identical observation likelihoods, although no pair of
    # individual class profiles is equal. No individual source is established.
    likelihood=np.array([[1,0,1,0],[0,1,0,1],[1,0,0,1],[0,1,1,0]],dtype=float)*np.exp(5.)
    assert np.array_equal(likelihood@np.array([.5,.5,0.,0.]),
        likelihood@np.array([0.,0.,.5,.5]))


def test_mixture_gradient_gap_matches_independent_dense_oracle():
    from vdjaqp.aggregate import fit_identity_mixture
    peptides=pl.DataFrame({'peptide_id':['p','q','r'],'spectrum_id':['s','t','u'],
        'seq':['ACDEFGHIK','YYYYYYYYY','RRRRRRRRR']})
    competition=pl.DataFrame({'spectrum_id':['s','s','t','u'],
        'identity_class':['A','B','B','A'],'log_weight':[3.,1.,4.,-1.]})
    out,receipt=fit_identity_mixture(competition,peptides,max_iterations=2)
    likelihood=np.column_stack([np.array([[np.exp(3.),np.exp(1.)],
        [0.,np.exp(4.)],[np.exp(-1.),0.]]),np.ones(3)])
    weights=np.r_[out['mixture_weight'].to_numpy(),receipt['null_weight']]
    denominator=likelihood@weights
    gradient=(likelihood/denominator[:,None]).sum(axis=0)
    assert np.isclose(receipt['max_component_gradient'],gradient.max())
    assert np.isclose(receipt['dual_gap'],max(0.,gradient.max()-3.))
    assert np.isclose(receipt['final_objective'],np.log(denominator).sum())
    # Every tested simplex point is bounded by the first-order concavity plane.
    upper=receipt['final_objective']+receipt['dual_gap']
    for a in np.linspace(0.,1.,21):
        for b in np.linspace(0.,1.-a,21):
            den=likelihood@np.array([a,b,1.-a-b])
            if (den>0).all():assert np.log(den).sum()<=upper+1e-12


def test_profiled_bounds_reject_four_class_alternative_explanations():
    from vdjaqp.aggregate import profile_identity_removals
    peptides=pl.DataFrame({'peptide_id':['p','q','r','t'],
        'spectrum_id':['V1','V2','J1','J2'],
        'seq':['ACDEFGHIK','YYYYYYYYY','RRRRRRRRR','PPPPPPPPP']})
    competition=pl.DataFrame({'spectrum_id':['V1','J1','V2','J2','V1','J2','V2','J1'],
        'identity_class':['A','A','B','B','C','C','D','D'],'log_weight':[5.]*8})
    out,receipt=profile_identity_removals(competition,peptides)
    assert receipt['performed_profiles']==4
    assert (out['fixed_other_removal_contrast']>.2).all()
    assert out['profile_score_available'].all()
    assert out['profiled_removal_lower_bound'].to_list()==[0.]*4
    assert any(not r['converged'] for r in receipt['reduced_fits'])


def test_profiled_bounds_retain_true_sources_and_explicit_partial_profiling():
    from vdjaqp.aggregate import profile_identity_removals
    peptides=pl.DataFrame({'peptide_id':['p','q'],'spectrum_id':['s','t'],
        'seq':['ACDEFGHIK','YYYYYYYYY']})
    competition=pl.DataFrame({'spectrum_id':['s','t'],'identity_class':['A','B'],
        'log_weight':[5.,3.]})
    out,receipt=profile_identity_removals(competition,peptides)
    assert (out['profiled_removal_lower_bound']>0.).all()
    assert (out['profiled_removal_lower_bound']<=out['profiled_removal_upper_bound']).all()
    partial,receipt=profile_identity_removals(competition,peptides,classes_to_profile=['A'])
    assert receipt['performed_profiles']==1
    assert partial.filter(pl.col('identity_class')=='B')['profile_status'][0]=='not_profiled'
    assert not partial.filter(pl.col('identity_class')=='B')['profile_score_available'][0]
    screened,receipt=profile_identity_removals(competition,peptides,minimum_contrast=100.)
    assert receipt['performed_profiles']==0
    assert screened['profiled_removal_lower_bound'].to_list()==[0.,0.]


def test_failed_reduced_numerical_certificate_is_unavailable(monkeypatch):
    import vdjaqp.aggregate as aggregate
    peptides=pl.DataFrame({'peptide_id':['p','q'],'spectrum_id':['s','t'],
        'seq':['ACDEFGHIK','YYYYYYYYY']})
    competition=pl.DataFrame({'spectrum_id':['s','t'],'identity_class':['A','B'],
        'log_weight':[5.,3.]})
    original=aggregate.fit_identity_mixture
    calls=0
    def fail_reduced(*args,**kwargs):
        nonlocal calls
        calls+=1
        if calls==1:return original(*args,**kwargs)
        raise ValueError('nonfinite mixture gradient/objective; no certified bound')
    monkeypatch.setattr(aggregate,'fit_identity_mixture',fail_reduced)
    out,receipt=aggregate.profile_identity_removals(competition,peptides)
    assert len(receipt['reduced_failures'])==2
    assert not out['profile_score_available'].any()
    assert out['profiled_removal_lower_bound'].null_count()==2
    assert out['profile_status'].to_list()==['reduced_fit_error']*2
