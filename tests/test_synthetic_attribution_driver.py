"""Constructed input/truth contracts only; no cohort search or annotation."""
import importlib.util
from pathlib import Path
import polars as pl
import pytest
from vdjaqp.readout import _SCHEMA
spec=importlib.util.spec_from_file_location('synthetic',Path(__file__).parents[1]/'scripts/compare_synthetic_attribution.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_sparse_zero_universe_and_true_providers():
    truth={('o1','a'),('o2','b')};predictions=[('o1','a',.2),('o1','a',.3),('o1','b',.1)]
    r=m.attribution_metrics(predictions,{'a','b','c'},{'o1','o2'},truth)
    assert r['n_pairs']==6 and r['true_positive']==1 and r['false_positive']==1 and r['false_negative']==1
    # Ranked groups: .5 correct, .1 wrong, zero tie of four with one positive.
    assert r['average_precision']==pytest.approx(.5+.5*2/6)
    assert m.attribution_metrics([],{'a'},{'light_only'},set())['average_precision'] is None
    with pytest.raises(ValueError):m.attribution_metrics([('o1','a',2)],{'a'},{'o1'},{('o1','a')})


def test_native_coisolation_and_light_only_truth(tmp_path):
    import json
    metadata=pl.DataFrame({'chain_id':['h1','h2','light'],'locus':['IGH','IGH','IGK'],
        'generator_junction_aa':['CAW','CLW','CFW']})
    peptides=pl.DataFrame({'enzyme':['Trypsin','Trypsin'],'native_spectrum_reference':['scan=1','scan=2'],
        'spectrum_id':['t:1','t:2']})
    path=tmp_path/'truth.json';path.write_text(json.dumps([
        dict(native_spectrum_reference='scan=1',contributors=[dict(peptide_aa='AAAAAA',protein_accessions=['h1','h2'])]),
        dict(native_spectrum_reference='scan=2',contributors=[dict(peptide_aa='AAAAAA',protein_accessions=['light'])])]))
    universe,observations,positive=m.native_truth_pairs(metadata,peptides,{'Trypsin':path})
    assert universe=={'CAW','CLW'} and observations=={'t:1','t:2'}
    assert positive=={('t:1','CAW'),('t:1','CLW')}


def test_input_projection_removes_old_labels(tmp_path):
    candidate={k:None for k in m.CANDIDATE_COLUMNS}
    candidate.update(clone_id='c',seq='CAAALFW',cdr3_seq='CAAALFW',cdr3_start=0,cdr3_end=7)
    pl.DataFrame([dict(candidate,background_member=False,confirmed_nats=999.,identity_status='secret_truth')]).write_parquet(tmp_path/'amplicon_candidate_identity.parquet')
    for enzyme in ('Trypsin','Chymotrypsin'):
        frame=pl.DataFrame(schema=_SCHEMA).with_columns(pl.lit(enzyme).alias('enzyme'),pl.lit(True).alias('correct'))
        frame.write_parquet(tmp_path/f'amplicon_{enzyme}_psm_correctness.parquet')
    import json,hashlib
    (tmp_path/'retained_table_manifest.json').write_text(json.dumps({'artifacts':[dict(file=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in tmp_path.glob('*.parquet')]}))
    candidates,peptides,_=m.load_inputs(tmp_path,'amplicon')
    assert candidates.height==1 and 'confirmed_nats' not in candidates.columns and 'identity_status' not in candidates.columns
    assert candidates.schema['background_member']==pl.Boolean and candidates['background_member'].to_list()==[False]
    assert 'correct' not in peptides.columns and 'native_contributor_count' not in peptides.columns


def test_one_native_annotation_batch_strips_probability_fields():
    candidates=pl.DataFrame({'clone_id':['c'],'seq':['CAAALFW'],'cdr3_start':[0],'cdr3_seq':['CAAALFW'],'v_name':['IGHV1'],'j_name':['IGHJ1']})
    calls=[]
    def annotate(a,v,j,**kw):
        calls.append((a,v,j,kw))
        return pl.DataFrame({'cdr3_repaired':['CAAALFW'],'cdr3_nt':['TGTGCCGCCGCCCTGTTTTGG'],
            'v_end_nt':[3],'j_start_nt':[18],'d_start_nt':[16],'d_end_nt':[18],
            'np1':['GCCGCCGCCCTG'],'np2':[None],'pgen':[1e-9]})
    result=m.annotate_candidates(candidates,2,annotator=annotate)
    assert len(calls)==1 and result.height==1 and 'pgen' not in result.columns


def test_j_frames_use_exact_scaffold_geometry_not_motifs():
    genes={'IGKJ1*01':'ACTGTT','IGLJ1*01':'GTTGTT','IGHJ1*01':'ACTGTT'}
    markup=pl.DataFrame({'scaffold_id':['k','l','h1','h2'],
        'j_call':['IGKJ1*01','IGLJ1*01','IGHJ1*01','IGHJ1*01'],
        'j_sequence_start':['2','3','1','2'],'vj_end':['7','8','6','7'],
        'fwr4_start':['4','4','1','4'],'fwr4_end':['7','8','6','7'],
        'fwr4':['TGTT','TTGTT','ACTGTT','TGTT']})
    frames,missing=m.j_scaffold_frames(genes,markup,{'k':'AACTGTT','l':'AAGTTGTT','h1':'ACTGTT','h2':'AACTGTT'})
    assert frames=={'IGKJ1*01':2,'IGLJ1*01':1}
    assert missing=={'IGHJ1*01':[0,2]}


def test_source_and_raw_pin_drift_is_rejected(tmp_path):
    path=tmp_path/'source.py';path.write_text('original')
    pins=m.file_hashes([path]);m.verify_pins(pins)
    path.write_text('changed')
    with pytest.raises(ValueError,match='changed'):m.verify_pins(pins)


def test_explicit_case_selection_is_unique_bounded_and_exclusive():
    assert m.selected_case_indices(0)==[0]
    assert m.selected_case_indices(case_indices=[1,3,17])==[1,3,17]
    assert m.selected_case_indices()==list(range(18))
    for index,indices in [(0,[1]),(None,[1,1]),(None,[-1]),(None,[18]),(None,[])]:
        with pytest.raises(ValueError):m.selected_case_indices(index,indices)


def test_operating_sweep_keeps_equal_scores_in_whole_groups():
    result,curve=m.attribution_metrics([('o1','a',.8),('o2','a',.4),('o1','b',.4),('o2','b',.2)],
        {'a','b'},{'o1','o2'},{('o1','a'),('o2','a')},return_curve=True)
    assert result['max_f1_threshold']==.4 and result['max_f1_f1']==pytest.approx(.8)
    point=[p for p in curve if p['threshold']==.4][0]
    assert (point['true_positive'],point['false_positive'],point['false_negative'])==(2,1,0)
    assert [p['n_called'] for p in curve]==[1,3,4]


def test_equal_max_f1_prefers_higher_threshold_and_fewer_calls():
    result,curve=m.attribution_metrics([('o1','a',.9),('o2','a',.5),('o1','b',.5),('o2','b',.5)],
        {'a','b'},{'o1','o2'},{('o1','a'),('o2','a')},return_curve=True)
    assert result['max_f1_threshold']==.9 and result['max_f1_n_called']==1
    assert curve[0]['f1']==curve[1]['f1']==pytest.approx(2/3)


def test_rounded_cutoff_is_used_for_counts_not_only_display():
    result,curve=m.attribution_metrics([('o','a',.67777777),('o','b',.66)],
        {'a','b'},{'o'},{('o','a')},return_curve=True)
    assert result['max_f1_threshold']==.67777777
    assert result['rounded_max_f1_threshold']==.67 and result['rounded_max_f1_f1']==1.
    fixed=m.operating_point(curve,.67,1)
    assert fixed['operating_true_positive']==1 and fixed['operating_false_positive']==0


def test_operating_parameters_require_training_provenance_and_are_fixed_on_test():
    import copy
    result,_=m.attribution_metrics([('o','a',.7)],{'a','b'},{'o'},{('o','a')},return_curve=True)
    rows=[dict(aggregation='alignment_union',view='unexplained_rescue',
        truth_definition='native_heavy_generative_N_full_codon_overlap',minimum_n_run=k,**result) for k in (0,2,3)]
    parameters=m.learn_operating_points(rows,{'source_sha256':{'fixture':'0'*64},'raw_input_sha256':{'fixture':'1'*64}})
    assert parameters['best_n_minimum']==2 and parameters['operating_points'][0]['threshold']==.7
    for mutate in [lambda d:d.update(training_case='test'),lambda d:d.update(raw_input_sha256={}),
                   lambda d:d['operating_points'][0].update(threshold=0),
                   lambda d:d['operating_points'].append(d['operating_points'][0])]:
        bad=copy.deepcopy(parameters);mutate(bad)
        with pytest.raises(ValueError):m.validate_operating_points(bad)
    _,curve=m.attribution_metrics([('o','a',.6),('o','b',.8)],{'a','b'},{'o'},{('o','a')},return_curve=True)
    fixed=m.operating_point(curve,parameters['operating_points'][0]['threshold'],1)
    assert fixed['operating_true_positive']==0 and fixed['operating_false_positive']==1


def test_mapped_positive_is_not_a_rescue_false_negative_and_k_universe_is_fixed():
    known={'mapped','unexplained','light_only'};positive={('mapped','a'),('unexplained','b')}
    all_predictions=[('mapped','a',.9),('unexplained','b',.6),('light_only','a',.2)]
    no_n=m.rescue_evaluation(all_predictions,{'a','b'},known,positive,{'mapped'})
    assert no_n['false_negative']==0 and no_n['n_positive_pairs']==1 and no_n['n_observations']==2
    for k_predictions in [[('unexplained','b',.6)],[]]:
        filtered=m.rescue_evaluation(k_predictions,{'a','b'},known,positive,{'mapped'})
        assert filtered['n_pairs']==no_n['n_pairs'] and filtered['n_positive_pairs']==no_n['n_positive_pairs']
        assert filtered['false_negative'] in (0,1)
    assert m.attribution_metrics(all_predictions,{'a','b'},known,positive)['n_positive_pairs']==2


def test_native_framework_and_coisolation_are_not_junction_or_n_truth(tmp_path):
    import json
    codons=dict(Q='CAA',C='TGT',A='GCC',S='TCT',L='CTG',F='TTT',W='TGG',E='GAA',G='GGC')
    physical=[];chains=[]
    for key,locus,variable in [('h1','IGH','QQCAAALFWEEGGGGGG'),('h2','IGH','QQCSSSLFWEEGGGGGG'),('light','IGK','QQCAAALFWEEGGGGGG')]:
        nt=''.join(codons[a] for a in variable);junction=variable[2:9]
        physical.append(dict(chain_id=key,locus=locus,variable_aa=variable,protein_aa=variable,
            generator_junction_aa=junction,generator_junction_nt=nt[6:27],generator_junction_start_nt=6,generator_junction_end_nt=27))
        a=dict(junction=nt[6:27],junction_aa=junction,junction_start=6,junction_end=27,
            v_sequence_end=9,d_sequence_start=18,d_sequence_end=24,j_sequence_start=24,
            np1=nt[9:18],np1_length=9,np1_n_tdt=9,np1_n_p=0,np1_n_inserted=0,
            np2='',np2_length=0,np2_n_tdt=0,np2_n_p=0,np2_n_inserted=0)
        chains.append(dict(chain_id=key,variable_nt=nt,variable_aa=variable,protein_aa=variable,airr_truth=a))
    metadata=pl.DataFrame(physical);generation={'chains':chains}
    peptides=pl.DataFrame({'enzyme':['Trypsin']*3,'native_spectrum_reference':['f','j','l'],'spectrum_id':['framework','junction','light_only']})
    path=tmp_path/'truth.json';path.write_text(json.dumps([
        dict(native_spectrum_reference='f',contributors=[dict(peptide_aa='EEGGGG',protein_accessions=['h1','h2'])]),
        dict(native_spectrum_reference='j',contributors=[dict(peptide_aa='QCAAAL',protein_accessions=['h1']),dict(peptide_aa='EEGGGG',protein_accessions=['h2'])]),
        dict(native_spectrum_reference='l',contributors=[dict(peptide_aa='EEGGGG',protein_accessions=['light'])])]))
    calls=[]
    def placements(heavy,queries,threads):
        calls.append((heavy.height,queries.height,threads));rows=[]
        for r in queries.iter_rows(named=True):
            for key in (['h1','h2'] if r['seq']=='EEGGGG' else ['h1']):
                start=9 if r['seq']=='EEGGGG' else 1
                rows.append(dict(chain_id=key,peptide_id=r['peptide_id'],start=start,end=start+6))
        return pl.DataFrame(rows)
    known,junction,n_known,n_positive,scope=m.native_junction_truth_pairs(metadata,peptides,{'Trypsin':path},generation,2,placement_fn=placements)
    assert calls==[(2,2,2)] and known==n_known=={'framework','junction','light_only'}
    assert junction==n_positive=={('junction','CAAALFW')}
    _,_,provider=m.native_truth_pairs(metadata,peptides,{'Trypsin':path})
    assert len(provider)==4 and ('framework','CSSSLFW') in provider
    assert scope['junction_unscorable_observations']==scope['strict_N_unscorable_observations']==0
    # A native provider peptide with literal copies inside and outside the junction
    # cannot identify which physical source location generated this spectrum.
    physical[0]['variable_aa']+='QCAAAL';physical[0]['protein_aa']+='QCAAAL'
    chains[0]['variable_aa']+='QCAAAL';chains[0]['protein_aa']+='QCAAAL'
    chains[0]['variable_nt']+=''.join(codons[a] for a in 'QCAAAL')
    def ambiguous_placements(heavy,queries,threads):
        extra=queries.filter(pl.col('seq')=='QCAAAL')['peptide_id'][0]
        return placements(heavy,queries,threads).vstack(pl.DataFrame([dict(chain_id='h1',peptide_id=extra,start=17,end=23)]))
    known,_,n_known,_,scope=m.native_junction_truth_pairs(pl.DataFrame(physical),peptides,{'Trypsin':path},generation,2,placement_fn=ambiguous_placements)
    assert known==n_known=={'framework','light_only'} and scope['junction_unscorable_observations']==1
    chains[0]['airr_truth'].update(np1_n_tdt=8,np1_n_p=1)
    _,_,n_known,n_positive,scope=m.native_junction_truth_pairs(pl.DataFrame(physical),peptides,{'Trypsin':path},generation,2,placement_fn=ambiguous_placements)
    assert n_known=={'light_only'} and not n_positive
    assert scope['strict_N_unscorable_observations']==2
