"""Evaluation-only metrics preserve coverage populations and contributor sets."""
import polars as pl
import pytest
from vdjaqp.simulation_metrics import exact_peptide_placements, region_coverage, psm_correctness


def test_exact_batch_repeats_and_short_rejection():
 chains=pl.DataFrame({'chain_id':['h'],'protein_aa':['ACDEFLACDEFI']})
 peptides=pl.DataFrame({'peptide_id':['p'],'seq':['ACDEFI']})
 hits=exact_peptide_placements(chains,peptides,threads=1)
 assert hits['start'].to_list()==[0,6]
 with pytest.raises(ValueError,match='unsearchable'):
  exact_peptide_placements(chains,peptides.with_columns(pl.lit('ACD').alias('seq')))


def test_coverage_union_zeros_and_unknown():
 regions=pl.DataFrame({'chain_id':['h','h','k','l'],'locus':['IGH','IGH','IGK','IGL'], 'region':['FR1','C','FR1','FR1'],'start':[0,5,0,None],'end':[5,10,3,None]})
 hits=pl.DataFrame({'chain_id':['h']*3,'start':[0,2,2],'end':[4,5,5],'evidence_level':['called']*3})
 out=region_coverage(regions,hits,['called','truth'])
 h=out.filter((pl.col('chain_id')=='h')&(pl.col('region')=='FR1')&(pl.col('evidence_level')=='called')).row(0,named=True)
 assert h['covered_residues']==5 and h['coverage']==1
 assert out.filter(pl.col('chain_id')=='k')['covered_residues'].to_list()==[0,0]
 assert out.filter(pl.col('chain_id')=='l')['coverage'].to_list()==[None,None]
 with pytest.raises(ValueError):region_coverage(regions,hits.with_columns(pl.lit(11).alias('end')),['called'])


def test_psm_sets_il_empty_and_missing():
 calls=pl.DataFrame({'peptide_id':['p','q','r','s'],'native_spectrum_reference':['a','a','b','missing'],'seq':['ACDEFI','YYYYYY','ACDEFL','ACDEFL']})
 truth=[{'native_spectrum_reference':'a','contributors':[{'peptide_aa':'ACDEFL','protein_accessions':['h','k']},{'peptide_aa':'GGGGGG','protein_accessions':['l']}]},{'native_spectrum_reference':'b','contributors':[{'peptide_aa':'ACDEFI','protein_accessions':['h']}]}]
 per,summary=psm_correctness(calls,truth)
 assert per['correct'].to_list()==[True,False,True,None]
 assert per['matching_protein_count'].to_list()==[2,0,1,None]
 assert summary.filter(pl.col('quantity')=='evaluated_psm_assignments')['count'].item()==3
 assert summary.filter(pl.col('quantity')=='evaluated_spectra')['count'].item()==2
 empty,stats=psm_correctness(calls.head(0),truth)
 assert empty.height==0 and stats.filter(pl.col('quantity')=='psm_assignments')['count'].item()==0
 with pytest.raises(ValueError):psm_correctness(calls,truth+[truth[0]])


def test_empty_coverage_and_input_keys():
 regions=pl.DataFrame({'chain_id':['h'],'locus':['IGH'],'region':['FR1'],'start':[0],'end':[4]})
 hits=pl.DataFrame(schema={'chain_id':pl.String,'start':pl.Int64,'end':pl.Int64,'evidence_level':pl.String})
 assert region_coverage(regions,hits,['called'])['covered_residues'].item()==0
 with pytest.raises(ValueError):region_coverage(pl.concat([regions,regions]),hits,['called'])
 calls=pl.DataFrame({'peptide_id':['p'],'native_spectrum_reference':[None],'seq':['ACDEFI']})
 with pytest.raises(ValueError,match='spectrum'):psm_correctness(calls,[])


def test_region_missing_reason_survives():
 regions=pl.DataFrame({'chain_id':['k'],'locus':['IGK'],'region':['D'],'start':[None],'end':[None],'status':['not_applicable'],'gene_call':[None]})
 hits=pl.DataFrame(schema={'chain_id':pl.String,'start':pl.Int64,'end':pl.Int64,'evidence_level':pl.String})
 out=region_coverage(regions,hits,['called'])
 assert out['region_status'].item()=='not_applicable'
 assert out['status'].item()=='undefined' and out['coverage'].item() is None


def test_identity_classes_full_partial_and_threshold():
 from vdjaqp.simulation_metrics import identity_metrics
 meta=pl.DataFrame({'chain_id':['h','h_alias','other','absent','k'],'locus':['IGH']*4+['IGK'],'variable_aa':['ACDEFI','ACDEFL','GGGGGG','YYYYYY','ACDEFI'],'variable_length':[6]*5,'protein_aa':['ACDEFIQQ','ACDEFLQQ','GGGGGGQQ','YYYYYYQQ','ACDEFIQQ']})
 candidates=pl.DataFrame({'clone_id':['full','alias','partial','none']})
 ranks=pl.DataFrame({'clone_id':candidates['clone_id'],'rank':[1,2,3,4],'n_tied':[1]*4,'confirmed_nats':[20.,21.,20.,19.9]})
 assoc=pl.DataFrame({'clone_id':['full','full','alias','partial','partial'],'chain_id':['h','h_alias','h','h','other'],'start':[0,0,0,1,1],'end':[6,6,6,4,4]})
 population,called,summary=identity_metrics(meta,candidates,ranks,assoc)
 assert population.height==4 and population['variable_mass_class_id'].n_unique()==3
 assert population.filter(pl.col('chain_id')=='h')['class_recovered_at20'].item()
 assert called.filter(pl.col('clone_id')=='partial')['identity_status'].item()=='partial_only'
 assert summary.filter(pl.col('quantity')=='recovered_classes_at20')['count'].item()==1
 assert population.filter(pl.col('chain_id')=='absent')['full_variable_compatible_candidate_count'].item()==0
 empty=assoc.head(0)
 _,_,s=identity_metrics(meta,candidates,ranks,empty)
 assert s.filter(pl.col('quantity')=='conditional_class_recovery_at20')['fraction'].item() is None
 with pytest.raises(ValueError):identity_metrics(meta,candidates,ranks.head(1),assoc)


def test_identity_ambiguous_full_and_empty_population():
 from vdjaqp.simulation_metrics import identity_metrics
 m=pl.DataFrame({'chain_id':['a','b'],'locus':['IGH','IGH'],'variable_aa':['ACDEFI','ACDEFIQ'],'variable_length':[6,7],'protein_aa':['ACDEFIQQ','ACDEFIQQ']})
 c=pl.DataFrame({'clone_id':['x']});r=pl.DataFrame({'clone_id':['x'],'rank':[1],'n_tied':[1],'confirmed_nats':[20.]})
 a=pl.DataFrame({'clone_id':['x','x'],'chain_id':['a','b'],'start':[0,0],'end':[8,8]})
 pop,called,stats=identity_metrics(m,c,r,a)
 assert called['identity_status'].item()=='ambiguous_full'
 assert not pop['class_recovered_at20'].any()
 pop,called,stats=identity_metrics(m.head(0),c.head(0),r.head(0),a.head(0))
 assert pop.height==called.height==0
 assert stats.filter(pl.col('quantity')=='recovered_classes_at20')['fraction'].item() is None
 with pytest.raises(ValueError):identity_metrics(m,c,r.with_columns(pl.lit(float('nan')).alias('confirmed_nats')),a)

 with pytest.raises(ValueError):identity_metrics(m,c.with_columns(pl.lit(5).alias('rank')),r,a)
