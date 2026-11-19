"""Optional emission components separate observed compatibility from context evidence."""
from types import SimpleNamespace
import numpy as np
import polars as pl
import pytest
from vdjaqp.alphabet import encode_frame
from vdjaqp.errormodel import ErrorModel,mass_prior
from vdjaqp.evidence import AlignmentScorer


def fixture(context_valid=True):
    candidates=pl.DataFrame({'clone_id':['exact','wrong'],'seq':['ACDEFG','ACDYFG']})
    code,offset,length=encode_frame(candidates['seq'])
    def frequency(slots,called):return np.where(np.asarray(slots)<6,.95,.001)
    background=SimpleNamespace(code=code,offset=offset,length=length,
        context_valid=np.full(code.size,context_valid,dtype=bool),
        frequency=frequency,mixed=lambda slots,called,matrix:frequency(slots,called))
    model=ErrorModel(np.array([0.,1.]),np.array([.1,.01]),mass_prior(),.1)
    scorer=AlignmentScorer(candidates,model,background)
    matches=pl.DataFrame({'clone_id':['wrong','exact'],'peptide_id':['p','p'],
        'start':[0,0],'end':[6,6],'n_subs':[1,0]})
    peptides=pl.DataFrame({'peptide_id':['p'],'spectrum_id':['s'],'seq':['ACDEFG'],
        'local_confidence':[[float('nan')]*6],'search_type':['db']})
    return scorer,matches,peptides


@pytest.mark.parametrize('chunk_size',[1,2,100])
def test_optional_components_preserve_default_llr_but_rank_observed_emission(monkeypatch,chunk_size):
    import vdjaqp.evidence as evidence
    monkeypatch.setattr(evidence,'_CHUNK',chunk_size)
    scorer,matches,peptides=fixture()
    default=scorer.score(matches,peptides)
    components=scorer.score(matches,peptides,include_components=True)
    assert components.select(default.columns).equals(default)
    assert np.allclose(components['llr'],components['log_context_emission']-components['log_context_null'])
    assert components['clone_id'].to_list()==matches['clone_id'].to_list()
    assert default['llr'][0]>default['llr'][1]
    assert components['log_emission'][0]<components['log_emission'][1]
    assert set(components.columns)-set(default.columns)=={'log_emission','log_context_emission','log_context_null'}


def test_unknown_context_does_not_erase_observed_mismatch():
    scorer,matches,peptides=fixture(False)
    components=scorer.score(matches,peptides,include_components=True)
    assert components['llr'].to_list()==[0.,0.]
    assert components['log_context_emission'].to_list()==[0.,0.]
    assert components['log_context_null'].to_list()==[0.,0.]
    assert components['log_emission'][0]<components['log_emission'][1]<0.
    assert np.isfinite(components['log_emission'].to_numpy()).all()


def test_optional_component_empty_schema():
    scorer,matches,peptides=fixture()
    result=scorer.score(matches.clear(),peptides,include_components=True)
    assert result.height==0
    assert all(result.schema[c]==pl.Float64 for c in
        ['llr','log_emission','log_context_emission','log_context_null'])


def test_evidence_from_frames_retains_components_in_single_build(monkeypatch):
    from vdjaqp import background as bg
    from vdjaqp.pipeline import Evidence
    monkeypatch.setattr(bg,'canonical_v',lambda:{})
    monkeypatch.setattr(bg,'canonical_j',lambda:{})
    clones=pl.DataFrame({'clone_id':['a','b'],'seq':['ACDEFG','ACDYFG'],
        'cdr3_start':[0,0],'cdr3_end':[6,6],'region_ids':[[6]*6,[6]*6],
        'v_name':['IGHV1-2']*2,'j_name':['IGHJ4']*2,'c_name':['IGHG1']*2})
    _,_,peptides=fixture()
    peptides=peptides.with_columns(pl.lit(.9).alias('alc'),pl.lit(1.).alias('alc_rank'),pl.lit(1).alias('n_reads'))
    model=ErrorModel(np.array([0.,1.]),np.array([.1,.01]),mass_prior(),.1)
    default=Evidence.from_frames(clones,peptides,error_model=model,max_subs=1)
    extended=Evidence.from_frames(clones,peptides,error_model=model,max_subs=1,include_components=True)
    assert extended.matches.equals(default.matches)
    assert extended.alignments.select(default.alignments.columns).equals(default.alignments)
    assert extended.score().equals(default.score())


def test_noncanonical_observed_residue_does_not_create_incomplete_emission():
    scorer,matches,peptides=fixture()
    scorer.background.code[0]=-1
    result=scorer.score(matches,peptides,include_components=True)
    assert np.isnan(result.filter(pl.col('clone_id')=='exact')['log_emission'][0])
    assert np.isfinite(result.filter(pl.col('clone_id')=='wrong')['log_emission'][0])


def test_common_emission_null_restores_compatibility_order_and_alias_invariance():
    from vdjaqp.aggregate import compete_identity_emissions
    scorer,matches,peptides=fixture()
    components=scorer.score(matches,peptides,include_components=True)
    metadata=peptides.with_columns(pl.lit(1.).alias('alc'),pl.lit(1.).alias('alc_rank'),pl.lit(1).alias('n_reads'))
    identities=pl.DataFrame({'clone_id':['wrong','exact','unmatched'],
        'identity_class':['B','A','C']})
    result,null=compete_identity_emissions(components,metadata,identities)
    exact=result.filter(pl.col('identity_class')=='A')['log_weight'][0]
    wrong=result.filter(pl.col('identity_class')=='B')['log_weight'][0]
    assert exact>wrong
    expected=np.log(np.exp(components['log_emission'].to_numpy()).sum()/3.)
    assert np.isclose(null['common_null_log_emission'][0],expected)
    alias=components.filter(pl.col('clone_id')=='exact').with_columns(pl.lit('alias').alias('clone_id'))
    other,other_null=compete_identity_emissions(pl.concat([components,alias]),metadata,
        pl.concat([identities,pl.DataFrame({'clone_id':['alias'],'identity_class':['A']})]))
    assert result.equals(other) and null.equals(other_null)


def test_common_null_zero_emissions_remain_null_only_and_unavailable_fails():
    from vdjaqp.aggregate import compete_identity_emissions
    _,_,peptides=fixture()
    metadata=peptides.with_columns(pl.lit(1.).alias('alc'),pl.lit(1.).alias('alc_rank'),pl.lit(1).alias('n_reads'))
    identities=pl.DataFrame({'clone_id':['a'],'identity_class':['A']})
    components=pl.DataFrame({'clone_id':['a'],'peptide_id':['p'],'log_emission':[float('-inf')]})
    result,null=compete_identity_emissions(components,metadata,identities)
    assert result.height==0
    assert null['supporting_classes'][0]==0
    assert null['null_status'][0]=='null_only_no_supported_class'
    assert null['common_null_log_emission'][0] is None
    with pytest.raises(ValueError,match='unavailable'):
        compete_identity_emissions(components.with_columns(pl.lit(float('nan')).alias('log_emission')),metadata,identities)


def test_emission_is_independent_of_candidate_context_background():
    scorer,matches,peptides=fixture()
    first=scorer.score(matches,peptides,include_components=True)
    scorer.background.frequency=lambda slots,called:np.full(np.asarray(slots).shape,.1)
    scorer.background.mixed=lambda slots,called,matrix:np.full(np.asarray(slots).shape,.1)
    second=scorer.score(matches,peptides,include_components=True)
    assert first['log_emission'].equals(second['log_emission'])
    assert not first['llr'].equals(second['llr'])


@pytest.mark.parametrize('log_emission',[0.,float('-inf')])
def test_unknown_observation_id_is_rejected_even_for_zero_emission(log_emission):
    from vdjaqp.aggregate import compete_identity_emissions
    _,_,peptides=fixture()
    metadata=peptides.with_columns(pl.lit(1.).alias('alc'),pl.lit(1.).alias('alc_rank'),pl.lit(1).alias('n_reads'))
    identities=pl.DataFrame({'clone_id':['a'],'identity_class':['A']})
    components=pl.DataFrame({'clone_id':['a'],'peptide_id':['unknown'],'log_emission':[log_emission]})
    with pytest.raises(ValueError,match='observation metadata'):
        compete_identity_emissions(components,metadata,identities)


@pytest.mark.parametrize('start,end', [(-1,5),(1,7),(0,5),(5,0),(None,6)])
def test_invalid_placement_cannot_read_another_reference(start,end):
    scorer,matches,peptides=fixture()
    malformed=matches.with_columns(pl.lit(start,dtype=pl.Int64).alias('start'),
                                   pl.lit(end).alias('end'))
    with pytest.raises(ValueError,match='placement span|reference bounds'):
        scorer.score(malformed,peptides)


def test_default_alignment_schema_stays_public():
    scorer,matches,peptides=fixture()
    assert scorer.score(matches,peptides).columns==matches.columns+['llr']
