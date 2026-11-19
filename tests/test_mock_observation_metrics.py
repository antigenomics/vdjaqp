"""Constructed evaluation-only checks; never invoke native annotation/search."""
import gzip
import hashlib
import polars as pl
import pytest
from scripts.export_mock_observation_metrics import queries_and_aliases, observability, bounded_coverage, raw_truth_digest


def test_duplicate_sequences_preserve_aliases_short_queries_excluded():
    entries=pl.DataFrame({'entry_id':['a','b','short'],'seq':['ACDEFG','ACDEFG','ACD']})
    candidates=pl.DataFrame({'clone_id':['c'],'seq':['ACDEFG']})
    q,a,d=queries_and_aliases(entries,candidates)
    assert q.height==1 and a.height==3 and d.height==1
    assert a['peptide_id'].n_unique()==1


def test_observability_partial_full_aliases_without_fake_scores():
    meta=pl.DataFrame({'chain_id':['h1','h2'],'locus':['IGH','IGH'],'variable_aa':['ACDEFG','ACDEFG'],'variable_length':[6,6]})
    hits=pl.DataFrame({'chain_id':['h1','h2','h1'],'peptide_id':['q','q','p'],'start':[0,0,1],'end':[6,6,5]})
    aliases=pl.DataFrame({'observed_id':['full','part'],'seq':['ACDEFG','CDEF'],'peptide_id':['q','p'],'evidence_level':['heavy_candidate']*2})
    pop,pairs=observability(meta,hits,aliases)
    assert pop.height==2 and pop['variable_mass_class_id'].n_unique()==1
    assert pop['variable_mass_class_id'][0]=='mass_'+hashlib.sha256(b'ACDEFG').hexdigest()
    assert pop['full_variable_compatible_candidate_count'].to_list()==[1,1]
    assert pairs.height==3 and pairs['full_variable_compatible'].sum()==2
    assert 'confirmed_nats' not in pop.columns and 'class_recovered_at20' not in pop.columns


def test_coverage_chunk_invariance_and_budget():
    from vdjaqp.simulation_metrics import region_coverage
    regions=pl.DataFrame({'chain_id':['a','b'],'locus':['IGH','IGK'],'region':['V','V'],'start':[0,0],'end':[8,8],'status':['defined']*2})
    placements=pl.DataFrame({'chain_id':['a','a','b'],'start':[0,2,1],'end':[4,6,3],'evidence_level':['arda_observed_exact_compatibility']*3})
    full=region_coverage(regions,placements,['arda_observed_exact_compatibility'])
    assert bounded_coverage(regions,placements,['arda_observed_exact_compatibility'],chains_per_chunk=1,max_slots=100)[0].equals(full)
    with pytest.raises(ValueError,match='expansion'):bounded_coverage(regions,placements,['arda_observed_exact_compatibility'],chains_per_chunk=1,max_slots=5)


def test_raw_truth_hash_stream(tmp_path):
    p=tmp_path/'truth.json.gz';p.write_bytes(gzip.compress(b'{"chains":[]}'))
    assert raw_truth_digest(p)==(hashlib.sha256(b'{"chains":[]}').hexdigest(),13)


@pytest.mark.parametrize('geometry',['amplicon','rnaseq'])
def test_export_observation_only_single_batches_and_receipt_tamper(tmp_path,monkeypatch,geometry):
    import json
    import vdjaqp.simulation_truth as truth_api
    import vdjaqp.simulation_metrics as metrics_api
    from scripts.export_mock_observation_metrics import main,digest
    branch=tmp_path/'branch';obs=branch/'observed';obs.mkdir(parents=True)
    tables={'catalog':pl.DataFrame({'source_observation_id':[geometry+'_observation'],'locus':['IGH'],'status':['candidate'],'background_member':[True]}),
            'candidates':pl.DataFrame({'clone_id':[geometry+'_c'],'seq':['ACDEFG'],'background_member':[True]}),
            'search_entries':pl.DataFrame({'entry_id':[geometry+'_c_0',geometry+'_c_1'],'seq':['ACDEFG','ACDEFG'],'background_member':[True,True]})}
    for name,frame in tables.items():frame.write_parquet(obs/(name+'.parquet'))
    native=branch/'native.tsv';native.write_text('native ARDA source fixture')
    receipt=obs/'receipt.json'
    record={'status':'complete','source_format':'arda_fastq_native','geometry':geometry,'counts':{n:f.height for n,f in tables.items()},'source_input_sha256':{str(native):digest(native)},'code_sha256':{},'artifact_sha256':{f.name:digest(f) for f in obs.iterdir()}}
    receipt.write_text(json.dumps(record))
    raw=json.dumps({'chains':[{'chain_id':'truth_only','locus':'IGH'}],'heavy_recombination_model':{'sha256':'model'}}).encode();packed=tmp_path/'truth.json.gz';packed.write_bytes(gzip.compress(raw))
    validation=tmp_path/'validation.json';validation.write_text(json.dumps({'status':'complete','observed_chains':1,'locus_counts':{'IGH':1},'model_sha256':'model','truth':{'raw_sha256':hashlib.sha256(raw).hexdigest(),'raw_bytes':len(raw),'gzip_sha256':digest(packed)}}))
    calls=[]
    meta=pl.DataFrame({'chain_id':['truth_only'],'locus':['IGH'],'variable_aa':['ACDEFG'],'variable_length':[6],'protein_aa':['ACDEFGH'],'v_call':['V*01'],'j_call':['J*01'],'constant_nt_gene':['C'],'canonical_v_gene':['IGHV1'],'canonical_v_match_status':['unique_gene'],'canonical_v_gene_count':[1]})
    regions=pl.DataFrame({'chain_id':['truth_only'],'locus':['IGH'],'region':['V'],'start':[0],'end':[7],'status':['defined'],'gene_call':['V*01']})
    def annotate(chains,threads):calls.append(('annotate',len(chains)));return meta,regions
    def placements(chains,queries,threads):
        calls.append(('placements',queries.height));assert queries['seq'].to_list()==['ACDEFG']
        return pl.DataFrame({'chain_id':['truth_only'],'peptide_id':queries['peptide_id'],'start':[0],'end':[6]})
    monkeypatch.setattr(truth_api,'annotate_truth',annotate);monkeypatch.setattr(metrics_api,'exact_peptide_placements',placements)
    output=tmp_path/'output';argv=['--truth-gzip',str(packed),'--validation-receipt',str(validation),'--observed-receipt',str(receipt),'--output',str(output),'--geometry',geometry,'--cpus','1']
    main(argv)
    assert calls==[('annotate',1),('placements',1)]
    assert pl.read_parquet(output/'gene_associated_capture.parquet')['supported_records'].to_list()==[1,1,1]
    assert pl.read_parquet(output/'heavy_identity_population.parquet')['full_variable_compatible_candidate_count'].item()==1
    assert not any('score' in f.name or 'recovery' in f.name for f in output.iterdir())
    assert json.loads((output/'provenance.json').read_text())['observability']['full_variable_observable_classes']==1
    with pytest.raises(FileExistsError):main(argv)
    (obs/'catalog.parquet').write_bytes(b'changed')
    with pytest.raises(ValueError,match='hash mismatch'):main([str(tmp_path/'bad') if x==str(output) else x for x in argv])
    assert len(calls)==2


def test_constant_only_does_not_count_partial_variable_observability():
    from scripts.export_mock_observation_metrics import observability
    m=pl.DataFrame({'chain_id':['h'],'locus':['IGH'],'variable_aa':['ACDEFG'],'variable_length':[6]})
    h=pl.DataFrame({'chain_id':['h'],'peptide_id':['q'],'start':[6],'end':[12]})
    a=pl.DataFrame({'observed_id':['c'],'seq':['HIKLMN'],'peptide_id':['q'],'evidence_level':['heavy_candidate']})
    p,associations=observability(m,h,a)
    assert p['partial_compatible_candidate_count'].item()==0
    assert associations['variable_overlap_residues'].item()==0


def test_empty_queries_and_no_hits_preserve_zero_observability():
    e=pl.DataFrame(schema={'entry_id':pl.String,'seq':pl.String})
    c=pl.DataFrame(schema={'clone_id':pl.String,'seq':pl.String})
    q,a,d=queries_and_aliases(e,c)
    assert q.height==a.height==d.height==0
    m=pl.DataFrame({'chain_id':['h'],'locus':['IGH'],'variable_aa':['ACDEFG'],'variable_length':[6]})
    h=pl.DataFrame(schema={'chain_id':pl.String,'peptide_id':pl.String,'start':pl.Int32,'end':pl.Int32})
    p,assoc=observability(m,h,a)
    assert p['full_variable_compatible_candidate_count'].item()==0 and assoc.height==0


def test_module_load_is_numerically_lazy():
    import subprocess,sys
    subprocess.run([sys.executable,'-c',"import runpy,sys;runpy.run_path('scripts/export_mock_observation_metrics.py',run_name='import_check');assert 'numpy' not in sys.modules and 'polars' not in sys.modules"],check=True,timeout=10)


def test_mass_class_namespace_matches_existing_il_collapse():
    from vdjaqp.match import ms_collapse
    m=pl.DataFrame({'chain_id':['h1','h2'],'locus':['IGH','IGH'],'variable_aa':['ACDEFI','ACDEFL'],'variable_length':[6,6]})
    hits=pl.DataFrame(schema={'chain_id':pl.String,'peptide_id':pl.String,'start':pl.Int32,'end':pl.Int32})
    aliases=pl.DataFrame(schema={'observed_id':pl.String,'seq':pl.String,'peptide_id':pl.String,'evidence_level':pl.String})
    population,_=observability(m,hits,aliases)
    expected=['mass_'+hashlib.sha256(seq.encode()).hexdigest() for seq in ms_collapse(m['variable_aa'])]
    assert population['variable_mass_class_id'].to_list()==expected
    assert len(set(expected))==1


def test_native_uncorrected_and_corrected_junction_aliases_remain_labelled():
    entries=pl.DataFrame({'entry_id':['read','corrected'],'seq':['ACDEFG','CAAFGW'],'sequence_origin':['arda_observed_read','arda_corrected_junction'],'coverage_kind':['observed_variable_fragment','junction_only'],'source_observation_id':['read_view','junction_view'],'observed_read_id':['native_read',None],'parent_clone_id':['parent','parent'],'native_clone_row':[0,0],'duplicate_count':[3,3]})
    candidates=pl.DataFrame(schema={'clone_id':pl.String,'seq':pl.String})
    queries,aliases,dropped=queries_and_aliases(entries,candidates)
    assert queries.height==2 and dropped.height==0
    assert set(aliases['sequence_origin'])=={'arda_observed_read','arda_corrected_junction'}
    assert aliases.filter(pl.col('sequence_origin')=='arda_corrected_junction')['coverage_kind'].item()=='junction_only'
    assert set(aliases['parent_clone_id'])=={'parent'}


def test_canonical_native_receipt_path_escape_and_remote_source_pins(tmp_path):
    import json
    from scripts.export_mock_observation_metrics import verify_observations,digest
    tables={'catalog':pl.DataFrame({'source_observation_id':['native1'],'background_member':[True]}),'candidates':pl.DataFrame({'clone_id':['native1'],'background_member':[True]}),'search_entries':pl.DataFrame({'entry_id':['native1'],'background_member':[True]})}
    for n,f in tables.items():f.write_parquet(tmp_path/(n+'.parquet'))
    receipt=tmp_path/'receipt.json';data=dict(status='complete',source_format='arda_fastq_native',geometry='amplicon',counts={n:f.height for n,f in tables.items()},artifact_sha256={n+'.parquet':digest(tmp_path/(n+'.parquet'))for n in tables},source_input_sha256={'/unmounted-native-test/raw.tsv':'recorded_remote_digest'})
    receipt.write_text(json.dumps(data));r,_,path=verify_observations(receipt,'amplicon')
    assert path==tmp_path.resolve() and r['source_verification']['unavailable_source_pins']==data['source_input_sha256']
    data['artifact_sha256']['../escape.parquet']='not-read';receipt.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='hash mismatch'):verify_observations(receipt,'amplicon')
    data['artifact_sha256'].pop('../escape.parquet');data['source_format']='retired_observation_route';receipt.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='native ARDA'):verify_observations(receipt,'amplicon')
