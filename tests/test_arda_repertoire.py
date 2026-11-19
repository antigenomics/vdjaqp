"""Constructed native TSV contracts; no read simulation or cohort processing."""
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess

import polars as pl
import pytest
from arda.refbuild.translate import translate

REPO=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('arda_repertoire',REPO/'scripts/smoke_repertoire.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def native_fixture(root, *, assembled=True):
    prefix=root/'native';nt='GCTTGTGCTTGGTGGGGTGGTTCTTCT';junction=nt[3:15]
    clones=pl.DataFrame({'locus':['IGH','IGK','IGL'],'junction':[junction]*3,
        'junction_aa':['CAWW']*3,'v_call':['IGHV1-2','IGKV1-5','IGLV1-40'],
        'j_call':['IGHJ4','IGKJ1','IGLJ1'],'c_call':['IGHG1','IGKC','IGLC1'],
        'duplicate_count':['7','2','3']})
    clones.write_csv(str(prefix)+'.clones.tsv',separator='\t')
    row={'sequence_id':'read0','sequence':nt,'rev_comp':'F','locus':'IGH',
         'junction':junction,'junction_aa':'CAWW','v_call':'IGHV1-2','j_call':'IGHJ4',
         'c_call':'IGHG1','fwr3_start':'1','fwr3_end':'6','cdr3_start':'7',
         'cdr3_end':'12','fwr4_start':'13','fwr4_end':'21','mmseqs2_score':'100'}
    # The unrelated row must be filtered by native clone keys before translation.
    rows=[row,{**row,'sequence_id':'read1','sequence':nt[:-3]+'GCT'},
          {**row,'sequence_id':'unrelated','v_call':'IGHV3-23','sequence':'N'*len(nt)}]
    pl.DataFrame(rows).write_csv(str(prefix)+'.airr.tsv',separator='\t')
    if assembled:
        pl.DataFrame([{key:value for key,value in {**row,'sequence_id':'rescued0',
                      'c_call':'','c_class':'IGHG'}.items()
                      if key in ['sequence_id','locus','v_call','j_call','c_call','c_class',
                                 'junction','junction_aa']}]).write_csv(str(prefix)+'.assembled.airr.tsv',separator='\t')
    pl.DataFrame({'sequence_id':['read0','read1'],'clone_row':[0,0]}).write_csv(str(prefix)+'.read_ids.tsv',separator='\t')
    return prefix


def test_native_selection_roles_and_coding_contract(tmp_path):
    prefix=native_fixture(tmp_path)
    report=module.export_candidates(prefix,tmp_path/'observed',geometry='rnaseq')
    frames={name:pl.read_parquet(tmp_path/'observed'/f'{name}.parquet')
            for name in ('catalog','candidates','search_entries')}
    assert report['counts']=={'catalog':3,'candidates':1,'search_entries':1}
    assert report['observed_locus_counts']=={'IGH':1,'IGK':1,'IGL':1}
    heavy=frames['candidates'].row(0,named=True)
    assert heavy['sequence_origin']=='arda_observed_read'
    assert heavy['observed_read_id']=='read0'
    assert heavy['native_matching_observations']==2
    assert heavy['native_sequence_variants']==2
    assert heavy['duplicate_count']==7 and heavy['native_clone_row']==0
    assert heavy['seq']=='ACAWWGG' and translate(heavy['nt_seq'])==heavy['seq']
    assert len(heavy['nt_seq'])==3*len(heavy['seq'])
    assert heavy['seq'][heavy['cdr3_start']:heavy['cdr3_end']]==heavy['native_junction_aa']
    assert frames['catalog']['source_observation_id'].n_unique()==3
    assert frames['catalog']['sequence_origin'].to_list()==[
        'arda_observed_read','arda_corrected_junction','arda_corrected_junction']
    for name,frame in frames.items():
        assert frame.schema['background_member']==pl.Boolean
        assert frame['background_member'].all()
        assert all(s.startswith('rnaseq_') for s in frame['source_observation_id'])
        assert report['artifact_sha256'][name+'.parquet']==module.digest(tmp_path/'observed'/f'{name}.parquet')
        assert 'source_contig_id' not in frame.columns
    assert len(report['source_input_sha256'])==4
    assert all(module.digest(p)==h for p,h in report['source_input_sha256'].items())
    with pytest.raises(FileExistsError):module.export_candidates(prefix,tmp_path/'observed')


@pytest.mark.parametrize('policy',['longest','native_root'])
def test_all_locus_candidates_preserve_heavy_outputs_and_observed_light_sequence(tmp_path,policy):
    prefix=native_fixture(tmp_path,assembled=False)
    reads=pl.read_csv(str(prefix)+'.airr.tsv',separator='\t',infer_schema_length=0)
    light=reads.head(1).with_columns(pl.lit('light').alias('sequence_id'),
        pl.lit('IGK').alias('locus'),pl.lit('IGKV1-5').alias('v_call'),
        pl.lit('IGKJ1').alias('j_call'),pl.lit('IGKC').alias('c_call'))
    pl.concat([reads,light]).write_csv(str(prefix)+'.airr.tsv',separator='\t')
    pl.DataFrame({'sequence_id':['read0','read1','light'],'clone_row':[0,0,1]}).write_csv(str(prefix)+'.read_ids.tsv',separator='\t')
    module.export_candidates(prefix,tmp_path/'default',representative_policy=policy)
    report=module.export_candidates(prefix,tmp_path/'all',all_loci=True,representative_policy=policy)
    assert not (tmp_path/'default/all_locus_candidates.parquet').exists()
    for name in ['catalog','candidates','search_entries']:
        assert (tmp_path/f'default/{name}.parquet').read_bytes()==(tmp_path/f'all/{name}.parquet').read_bytes()
    candidates=pl.read_parquet(tmp_path/'all/all_locus_candidates.parquet')
    assert candidates['locus'].to_list()==['IGH','IGK','IGL']
    assert candidates['clone_id'].n_unique()==3
    heavy=pl.read_parquet(tmp_path/'all/candidates.parquet')
    assert candidates.filter(pl.col('locus')=='IGH').equals(heavy)
    kappa=candidates.filter(pl.col('locus')=='IGK').row(0,named=True)
    assert kappa['seq']=='ACAWWGG' and translate(kappa['nt_seq'])==kappa['seq']
    assert kappa['v_name']=='IGKV1-5' and kappa['observed_read_id']=='light'
    assert len(kappa['region_ids'])==len(kappa['seq'])
    lam=candidates.filter(pl.col('locus')=='IGL').row(0,named=True)
    assert lam['coverage_kind']=='junction_only' and lam['seq']=='CAWW'
    assert lam['coding_nt_start'] is None and lam['observed_read_id'] is None
    assert report['counts']['all_locus_candidates']==3
    assert report['artifact_sha256']['all_locus_candidates.parquet']==module.digest(tmp_path/'all/all_locus_candidates.parquet')


def test_late_nullable_observed_call_keeps_string_schema(tmp_path):
    """A valid C call after 100 junction-only roots must not infer a Null column."""
    prefix=native_fixture(tmp_path,assembled=False)
    heavy=pl.read_csv(str(prefix)+'.clones.tsv',separator='\t').row(0,named=True)
    roots=[{**heavy,'v_call':f'IGHV1-{i+10}'} for i in range(101)]+[heavy]
    pl.DataFrame(roots).write_csv(str(prefix)+'.clones.tsv',separator='\t')
    pl.DataFrame({'sequence_id':['read0','read1'],'clone_row':[101,101]}).write_csv(str(prefix)+'.read_ids.tsv',separator='\t')
    module.export_candidates(prefix,tmp_path/'observed')
    for name in ('catalog','candidates','search_entries'):
        frame=pl.read_parquet(tmp_path/'observed'/f'{name}.parquet')
        assert frame.schema['observed_c_call']==pl.String
        assert frame.filter(pl.col('observed_c_call').is_not_null())['observed_c_call'].to_list()==['IGHG1']


def test_explicit_gzip_airr_preserves_complete_exports_and_actual_source_pins(tmp_path):
    prefix=native_fixture(tmp_path)
    raw=Path(str(prefix)+'.airr.tsv');packed=raw.with_suffix('.tsv.gz')
    packed.write_bytes(gzip.compress(raw.read_bytes(),mtime=0))
    baseline=module.export_candidates(prefix,tmp_path/'raw_observed')
    fresh=module.export_candidates(prefix,tmp_path/'gzip_observed',raw_airr=packed)
    for name in ['catalog','candidates','search_entries']:
        assert pl.read_parquet(tmp_path/'raw_observed'/f'{name}.parquet').equals(
            pl.read_parquet(tmp_path/'gzip_observed'/f'{name}.parquet'))
    assert baseline['counts']==fresh['counts']
    assert str(packed.resolve()) in fresh['source_input_sha256']
    assert str(raw.resolve()) not in fresh['source_input_sha256']
    assert fresh['source_input_sha256'][str(packed.resolve())]==module.digest(packed)
    raw.unlink()
    with pytest.raises(FileNotFoundError):module.export_candidates(prefix,tmp_path/'implicit_missing')
    assert module.export_candidates(prefix,tmp_path/'gzip_only',raw_airr=packed)['counts']==fresh['counts']
    with pytest.raises(FileNotFoundError):module.export_candidates(prefix,tmp_path/'explicit_missing',raw_airr=tmp_path/'absent.tsv.gz')


def test_all_loci_observed_entries_and_longest_frame(tmp_path):
    prefix=native_fixture(tmp_path,assembled=False)
    raw=pl.read_csv(str(prefix)+'.airr.tsv',separator='\t')
    row=raw.row(0,named=True)
    extra=[{**row,'sequence_id':'short','sequence':row['sequence'][:18], 'fwr4_end':18},
           {**row,'sequence_id':'kappa','locus':'IGK','v_call':'IGKV1-5','j_call':'IGKJ1'},
           {**row,'sequence_id':'lambda','locus':'IGL','v_call':'IGLV1-40','j_call':'IGLJ1'}]
    pl.concat([raw,pl.DataFrame(extra)],how='vertical_relaxed').write_csv(str(prefix)+'.airr.tsv',separator='\t')
    pl.DataFrame({'sequence_id':['read0','read1','short','kappa','lambda'],'clone_row':[0,0,0,1,2]}).write_csv(str(prefix)+'.read_ids.tsv',separator='\t')
    report=module.export_candidates(prefix,tmp_path/'observed')
    entries=pl.read_parquet(tmp_path/'observed/search_entries.parquet')
    candidate=pl.read_parquet(tmp_path/'observed/candidates.parquet').row(0,named=True)
    assert entries['locus'].to_list()==['IGH','IGK','IGL']
    assert report['counts']=={'catalog':3,'candidates':1,'search_entries':3}
    assert candidate['observed_read_id']=='read0' and candidate['seq']=='ACAWWGG'
    assert candidate['native_matching_observations']==3
    assert candidate['sequence_origin']=='arda_observed_read'


def test_shorter_assembled_compatibility_frame_cannot_replace_longer_observation(tmp_path):
    prefix=native_fixture(tmp_path,assembled=False)
    row=pl.read_csv(str(prefix)+'.airr.tsv',separator='\t').row(0,named=True)
    # Compatibility with a version supplying full sequences is retained, but
    # current native Stage3 emits junction fields only (the primary fixture).
    pl.DataFrame([{**row,'sequence_id':'assembled_short','sequence':row['sequence'][:18],
        'fwr4_end':18,'mmseqs2_score':10000}]).write_csv(
            str(prefix)+'.assembled.airr.tsv',separator='\t')
    pl.DataFrame({'sequence_id':['read0','read1','assembled_short'],'clone_row':[0,0,0]}).write_csv(str(prefix)+'.read_ids.tsv',separator='\t')
    report=module.export_candidates(prefix,tmp_path/'observed')
    candidate=pl.read_parquet(tmp_path/'observed/candidates.parquet').row(0,named=True)
    assert candidate['observed_read_id']=='read0'
    assert candidate['sequence_origin']=='arda_observed_read'
    assert candidate['seq']=='ACAWWGG' and len(candidate['nt_seq'])==21
    assert candidate['native_matching_observations']==3
    assert report['observed_variable_fragment_candidates']==1


@pytest.mark.parametrize('timeout',['0','-1','nan','inf'])
def test_cli_rejects_unbounded_timeout_before_outputs(tmp_path,timeout):
    with pytest.raises(SystemExit):module.main(['amplicon','--root',str(tmp_path),'--stage-timeout',timeout])
    assert not list(tmp_path.iterdir())


def test_preprocess_timeout_is_passed_without_filter_changes(tmp_path,monkeypatch):
    from scripts import smoke_preprocess
    calls=[]
    def fail(command,**kwargs):
        calls.append((command,kwargs));raise subprocess.TimeoutExpired(command,kwargs['timeout'])
    monkeypatch.setattr(smoke_preprocess.subprocess,'run',fail)
    with pytest.raises(subprocess.TimeoutExpired):
        smoke_preprocess.preprocess_pairs('r1','r2',tmp_path/'fresh.gz',threads=16,fastp='fastp',stage_timeout=12)
    command,kwargs=calls[0]
    assert kwargs=={'check':True,'timeout':12}
    assert '--merge' in command and '--include_unmerged' in command
    assert all(flag in command for flag in ['--disable_adapter_trimming','--disable_trim_poly_g',
        '--disable_quality_filtering','--disable_length_filtering'])


def test_cli_one_native_batch_default_modes_and_complete_receipt(tmp_path,monkeypatch):
    from scripts import smoke_preprocess
    import vdjaqp
    import shutil
    prefix=native_fixture(tmp_path)
    for mate in ['R1','R2']:
        with gzip.open(tmp_path/f'amplicon_{mate}.fastq.gz','wt') as stream:
            stream.write('@software_string\nACGT\n+\nIIII\n')
    budgets=[];native=[];preprocess=[]
    monkeypatch.setattr(vdjaqp,'resolve_cpu_budget',lambda requested:budgets.append(requested) or 16)
    def merge(*args,**kwargs):
        preprocess.append(kwargs);Path(args[2]).write_bytes(b'constructed');return {'input_pairs':1}
    def run(command,**kwargs):
        native.append((command,kwargs))
        out=Path(command[command.index('--out-dir')+1])
        for suffix in ['.clones.tsv','.airr.tsv','.assembled.airr.tsv']:
            shutil.copyfile(str(prefix)+suffix,str(out/'smoke')+suffix)
        pl.DataFrame({'sequence_id':['read0','read1'],'clone_row':[0,0]}).write_csv(str(out/'smoke')+'.read_ids.tsv',separator='\t')
    monkeypatch.setattr(smoke_preprocess,'preprocess_pairs',merge)
    monkeypatch.setattr(module.subprocess,'run',run)
    module.main(['amplicon','--root',str(tmp_path),'--threads','100','--stage-timeout','42'])
    assert budgets==[100] and len(native)==len(preprocess)==1
    command,kwargs=native[0]
    assert command[1]=='amplicon' and command[-2:]==['--threads','16']
    assert '--reconstruct' not in command and '--ec-mode' not in command and '--r2' not in command
    assert kwargs=={'check':True,'timeout':42.0} and preprocess[0]['stage_timeout']==42
    receipt=json.loads((tmp_path/'amplicon_observed/receipt.json').read_text())
    assert receipt['status']=='complete' and receipt['cpu_budget']==16
    assert len(receipt['input_fastq_sha256'])==2 and receipt['preprocess']['input_pairs']==1


def test_raw_clone_geometry_preserved_without_changing_candidate_or_coverage(tmp_path):
    prefix=native_fixture(tmp_path)
    module.export_candidates(prefix,tmp_path/'before')
    fields=dict(d_call='IGHD1-1*01,IGHD1-1*02',d2_call='IGHD2-2*01',
        d_support='1.20e-04',d2_support='2.30e-03',np1='A',np2='TG',np3='C',
        v_sequence_end='3',d_sequence_start='5',d_sequence_end='6',
        d2_sequence_start='9',d2_sequence_end='10',j_sequence_start='12')
    clones=pl.read_csv(str(prefix)+'.clones.tsv',separator='\t',infer_schema_length=0)
    clones.with_columns(*[pl.when(pl.col('locus')=='IGH').then(pl.lit(value)).otherwise(pl.lit('-1') if 'sequence_' in key else None).alias(key) for key,value in fields.items()]).write_csv(str(prefix)+'.clones.tsv',separator='\t')
    module.export_candidates(prefix,tmp_path/'after')
    for name in ('catalog','candidates','search_entries'):
        before=pl.read_parquet(tmp_path/'before'/f'{name}.parquet')
        after=pl.read_parquet(tmp_path/'after'/f'{name}.parquet')
        raw=[c for c in after.columns if c.startswith('native_clone_') and c!='native_clone_row']
        assert after.drop(raw).equals(before.drop(raw))
        heavy=after.filter(pl.col('locus')=='IGH').row(0,named=True)
        assert all(heavy['native_clone_'+key]==value for key,value in fields.items())
        assert all(after.schema['native_clone_'+key]==pl.String for key in fields)
        assert 'corrected-junction' in heavy['native_clone_geometry_provenance']
        assert not set(fields)&set(after.columns)
    catalog=pl.read_parquet(tmp_path/'after/catalog.parquet')
    assert catalog.filter(pl.col('locus')=='IGK')['native_clone_d_sequence_start'].item()=='-1'
    absent=pl.read_parquet(tmp_path/'before/candidates.parquet')
    assert all(absent['native_clone_'+key].null_count()==absent.height for key in fields)


@pytest.mark.parametrize('updates,expected,status,bounds',[
    ({'v_sequence_end':'4','j_sequence_start':'9'},[4,6,6,5,7],'inferred_vj_excluded_core',(4,8)),
    ({'v_sequence_end':'6','j_sequence_start':'7'},[4,5,5,5,7],'inferred_vj_excluded_core',(6,6)),
    ({},[4,6,6,6,7],'conservative_whole_junction_unresolved',(None,None)),
    ({'v_sequence_end':'-1','j_sequence_start':'10'},[4,6,6,6,7],'conservative_whole_junction_unresolved',(None,None)),
    ({'v_sequence_end':'3.0','j_sequence_start':'10'},[4,6,6,6,7],'conservative_whole_junction_unresolved',(None,None)),
    ({'v_sequence_end':'3','j_sequence_start':'14'},[4,6,6,6,7],'conservative_whole_junction_unresolved',(None,None)),
    ({'v_sequence_end':'3','j_sequence_start':'10','v_call':'IGHV1-2,IGHV3-23'},[4,6,6,6,7],'conservative_whole_junction_unresolved',(None,None)),
    ({'v_sequence_end':'3','j_sequence_start':'10','j_call':'unknown'},[4,6,6,6,7],'conservative_whole_junction_unresolved',(None,None)),
])
def test_junction_core_protection_geometry(updates,expected,status,bounds):
    row=dict(locus='IGH',v_call='IGHV1-2*01',j_call='IGHJ4*02',junction='TGTGCTGCTGCTTGG',junction_aa='CAAAW')
    row.update(updates);labels,metadata=module.junction_protection(row)
    assert labels==expected and labels[0]==4 and labels[-1]==7
    assert metadata['junction_core_status']==status
    assert (metadata['junction_core_start_nt'],metadata['junction_core_end_nt'])==bounds


def test_export_protection_survives_actual_shm_consumer(tmp_path,monkeypatch):
    from vdjaqp import shm
    prefix=native_fixture(tmp_path)
    module.export_candidates(prefix,tmp_path/'observed')
    candidates=pl.read_parquet(tmp_path/'observed/candidates.parquet')
    row=candidates.row(0,named=True);start,end=row['cdr3_start'],row['cdr3_end']
    assert row['region_ids'][start:end]==[4,6,6,7]
    assert row['junction_core_status']=='conservative_whole_junction_unresolved'
    monkeypatch.setattr(shm,'canonical_v',lambda:{'IGHV1-2':'A'*30})
    monkeypatch.setattr(shm,'canonical_j',lambda:{'IGHJ4':'A'*30})
    monkeypatch.setattr(shm,'aligned_columns',lambda query,reference:((i,i,'A')for i in range(len(query))))
    reverted=shm.revert_junctions(candidates)['junction_reverted'].item()
    assert reverted=='AAWA' and reverted[1:-1]==row['cdr3_seq'][1:-1]


def test_corrected_membership_retains_long_synonymous_observation_without_pasting(tmp_path):
    prefix=tmp_path/'native';native='TGTGCTTGG';raw='TGCGCCTGG'
    pl.DataFrame(dict(locus=['IGH'],junction=[native],junction_aa=['CAW'],v_call=['IGHV1-2'],
        j_call=['IGHJ4'],c_call=['IGHG1'],duplicate_count=['2'])).write_csv(str(prefix)+'.clones.tsv',separator='\t')
    common=dict(rev_comp='F',locus='IGH',v_call='IGHV1-2',j_call='IGHJ4',c_call='IGHG1',junction_aa='CAW',mmseqs2_score='100')
    sequence='GCT'*5+raw+'GGT'*8
    rows=[dict(common,sequence_id='long',sequence=sequence,junction=raw,v_call='IGHV3-23',fwr3_start='1',fwr3_end='18',cdr3_start='19',cdr3_end='21',fwr4_start='22',fwr4_end='48'),
          dict(common,sequence_id='short',sequence=native,junction=native,fwr3_start='1',fwr3_end='3',cdr3_start='4',cdr3_end='6',fwr4_start='7',fwr4_end='9')]
    pl.DataFrame(rows).write_csv(str(prefix)+'.airr.tsv',separator='\t')
    with pytest.raises(FileNotFoundError,match='membership is required'):
        module.export_candidates(prefix,tmp_path/'missing_map')
    assert not (tmp_path/'missing_map').exists()
    pl.DataFrame({'sequence_id':['long','short'],'locus':['IGH','IGH'],'junction':[native,native],'clone_row':[0,0]}).write_csv(str(prefix)+'.read_ids.tsv',separator='\t')
    module.export_candidates(prefix,tmp_path/'fixed')
    chosen=pl.read_parquet(tmp_path/'fixed/candidates.parquet').row(0,named=True)
    assert chosen['seq']=='AAAAACAWGGGGGGGG' and chosen['nt_seq']==sequence
    assert chosen['observed_junction_nt']==raw and chosen['native_junction_nt']==native
    assert chosen['cdr3_seq']==chosen['observed_junction_aa']=='CAW'
    assert chosen['native_clone_row']==0 and chosen['duplicate_count']==2
    assert chosen['native_matching_observations']==2
    assert chosen['observed_v_call']=='IGHV3-23' and chosen['v_name']=='IGHV1-2'
    assert chosen['native_membership_source']=='native_correct_airr_read_map'
    assert chosen['junction_core_status']=='observed_junction_differs_from_corrected_root'


def test_observed_junction_uses_native_nt_coordinate_when_peptide_repeats(tmp_path):
    junction='TGTGCTTGG'
    sequence='GCT'+'TGCGCCTGG'+'GGT'+junction+'GGTGGT'
    assert sequence.count(junction)==1 and translate(sequence).count('CAW')==2
    path=tmp_path/'observed.tsv'
    pl.DataFrame([dict(sequence_id='observed',sequence=sequence,rev_comp='F',
        locus='IGH',v_call='IGHV1-2',j_call='IGHJ4',c_call='IGHG1',
        junction=junction,junction_aa='CAW',mmseqs2_score='100',
        fwr3_start='1',fwr3_end='18',cdr3_start='19',cdr3_end='21',
        fwr4_start='22',fwr4_end='30')]).write_csv(path,separator='\t')
    projected=module._observations(path,False,pl.DataFrame({'sequence_id':['observed']})).collect()
    assert projected.height==1
    row=projected.row(0,named=True)
    assert row['seq']=='ACAWGCAWGG' and row['_js']==5
    assert row['_aa'][row['_js']:row['_js']+3]=='CAW'


@pytest.mark.parametrize('aa,nt',[('AAW','GCTGCTTGG'),('CAA','TGTGCTGCT'),('C*W','TGTTAGTGG'),('CAW','TGTGCTTG')])
def test_junction_protection_requires_complete_native_anchor_geometry(aa,nt):
    with pytest.raises(ValueError,match='complete'):
        module.junction_protection(dict(locus='IGH',v_call='IGHV1-2',j_call='IGHJ4',junction=nt,junction_aa=aa))


def test_partial_observation_is_retained_but_not_selected_as_full_junction_query(tmp_path):
    path=tmp_path/'partial.tsv'
    pl.DataFrame([dict(sequence_id='partial',sequence='GCTGCTTGGGGT',rev_comp='F',locus='IGH',v_call='IGHV1-2',j_call='IGHJ4',c_call='IGHG1',junction='GCTGCTTGG',junction_aa='AAW',fwr3_start='1',fwr3_end='3',cdr3_start='4',cdr3_end='6',fwr4_start='7',fwr4_end='12',mmseqs2_score='100')]).write_csv(path,separator='\t')
    frame=module._observations(path,False,pl.DataFrame({'sequence_id':['partial']})).collect()
    assert frame.height==1 and frame['native_complete_junction'].item() is False


@pytest.mark.parametrize('complete_observation',[True,False])
def test_native_export_prefers_complete_frame_or_declares_junction_fallback(tmp_path,complete_observation):
    prefix=native_fixture(tmp_path,assembled=False)
    path=Path(str(prefix)+'.airr.tsv')
    original=pl.read_csv(path,separator='\t',infer_schema_length=0)
    nt='GCT'*10+'GCTGCTTGG'+'GGT'*20
    partial=dict(sequence_id='long_partial',sequence=nt,rev_comp='F',locus='IGH',
        v_call='IGHV1-2',j_call='IGHJ4',c_call='IGHG1',junction='GCTGCTTGG',
        junction_aa='AAW',fwr3_start='1',fwr3_end='33',cdr3_start='34',
        cdr3_end='36',fwr4_start='37',fwr4_end=str(len(nt)),mmseqs2_score='200')
    frame=pl.concat([original,pl.DataFrame([partial])],how='diagonal_relaxed') if complete_observation else pl.DataFrame([partial])
    frame.write_csv(path,separator='\t')
    membership=pl.read_csv(str(prefix)+'.read_ids.tsv',separator='\t')
    link=pl.DataFrame({'sequence_id':['long_partial'],'clone_row':[0]})
    pl.concat([membership,link]).write_csv(str(prefix)+'.read_ids.tsv',separator='\t') if complete_observation else link.write_csv(str(prefix)+'.read_ids.tsv',separator='\t')
    module.export_candidates(prefix,tmp_path/'result')
    row=pl.read_parquet(tmp_path/'result/candidates.parquet').row(0,named=True)
    assert row['native_matching_observations']==(3 if complete_observation else 1)
    assert row['native_complete_matching_observations']==(2 if complete_observation else 0)
    if complete_observation:
        assert row['observed_read_id']=='read0' and row['sequence_origin']=='arda_observed_read'
    else:
        assert row['seq']=='CAWW' and row['sequence_origin']=='arda_corrected_junction'


def test_incomplete_junction_does_not_supply_missing_region_geometry(tmp_path):
    path=tmp_path/'unresolved.tsv'
    pl.DataFrame([dict(sequence_id='partial',sequence='GCTGCTTGGTGG',rev_comp='F',locus='IGH',v_call='IGHV1-2',j_call='IGHJ4',c_call='IGHG1',junction='GCTGCTTGG',junction_aa='AAW',fwr3_start='1',fwr3_end='3',fwr4_start='7',fwr4_end='12',mmseqs2_score='100')]).write_csv(path,separator='\t')
    assert module._observations(path,False,pl.DataFrame({'sequence_id':['partial']})).collect().height==0


def root_policy_fixture(root, *, include_root=True):
    prefix=native_fixture(root,assembled=False)
    path=Path(str(prefix)+'.airr.tsv')
    original=pl.read_csv(path,separator='\t',infer_schema_length=0)
    wrong=dict(sequence_id='long_wrong',sequence='GCT'*3+'TGTGTTTGGTGG'+'GGT'*4,
        rev_comp='F',locus='IGH',v_call='IGHV1-2',j_call='IGHJ4',c_call='IGHG1',
        junction='TGTGTTTGGTGG',junction_aa='CVWW',fwr3_start='1',fwr3_end='12',
        cdr3_start='13',cdr3_end='18',fwr4_start='19',fwr4_end='33',mmseqs2_score='500')
    rows=pl.concat([original,pl.DataFrame([wrong])],how='diagonal_relaxed') if include_root else pl.DataFrame([wrong])
    rows.write_csv(path,separator='\t')
    links=pl.DataFrame({'sequence_id':['long_wrong'],'clone_row':[0]})
    if include_root:links=pl.concat([pl.read_csv(str(prefix)+'.read_ids.tsv',separator='\t'),links])
    links.write_csv(str(prefix)+'.read_ids.tsv',separator='\t')
    return prefix


def test_native_root_policy_prefers_observed_matching_junction_without_population_filter(tmp_path):
    prefix=root_policy_fixture(tmp_path)
    default=module.export_candidates(prefix,tmp_path/'default')
    explicit=module.export_candidates(prefix,tmp_path/'explicit',representative_policy='longest')
    root=module.export_candidates(prefix,tmp_path/'root',representative_policy='native_root')
    assert default['counts']==explicit['counts']==root['counts']
    for name in ['catalog','candidates','search_entries']:
        before=pl.read_parquet(tmp_path/'default'/f'{name}.parquet')
        same=pl.read_parquet(tmp_path/'explicit'/f'{name}.parquet')
        after=pl.read_parquet(tmp_path/'root'/f'{name}.parquet')
        assert before.equals(same)
        assert before['background_member'].equals(after['background_member'])
        assert before['parent_clone_id'].to_list()==after['parent_clone_id'].to_list()
    old=pl.read_parquet(tmp_path/'default/candidates.parquet').row(0,named=True)
    chosen=pl.read_parquet(tmp_path/'root/candidates.parquet').row(0,named=True)
    assert old['observed_read_id']=='long_wrong' and old['cdr3_seq']=='CVWW'
    assert chosen['observed_read_id']=='read0' and chosen['cdr3_seq']=='CAWW'
    assert chosen['observed_junction_nt']==chosen['native_junction_nt']
    assert len(chosen['seq'])<len(old['seq'])
    assert chosen['representative_selection_status']=='observed_native_root'
    assert chosen['representative_native_root_match'] is True
    assert chosen['native_root_matching_observations']==2
    assert chosen['native_matching_observations']==old['native_matching_observations']==3
    assert translate(chosen['nt_seq'])==chosen['seq']
    assert 'representative_selection_status' not in old


def test_native_root_policy_falls_back_explicitly_without_inventing_sequence(tmp_path):
    prefix=root_policy_fixture(tmp_path,include_root=False)
    module.export_candidates(prefix,tmp_path/'default')
    report=module.export_candidates(prefix,tmp_path/'root',representative_policy='native_root')
    old=pl.read_parquet(tmp_path/'default/candidates.parquet')
    chosen=pl.read_parquet(tmp_path/'root/candidates.parquet')
    assert chosen.select(old.columns).equals(old)
    row=chosen.row(0,named=True)
    assert row['cdr3_seq']=='CVWW' and row['native_junction_aa']=='CAWW'
    assert row['representative_native_root_match'] is False
    assert row['representative_selection_status']=='observed_fallback_no_root_match'
    assert row['native_root_matching_observations']==0
    catalog=pl.read_parquet(tmp_path/'root/catalog.parquet')
    assert catalog.filter(pl.col('coverage_kind')=='junction_only')['representative_native_root_match'].null_count()==2
    assert report['representative_policy']=='native_root'


def test_native_root_policy_retains_deterministic_read_id_ties(tmp_path):
    prefix=native_fixture(tmp_path,assembled=False)
    row=pl.read_csv(str(prefix)+'.airr.tsv',separator='\t',infer_schema_length=0).row(0,named=True)
    for folder,ids in [('one',['zeta','alpha']),('two',['alpha','zeta'])]:
        pl.DataFrame([{**row,'sequence_id':name} for name in ids]).write_csv(str(prefix)+'.airr.tsv',separator='\t')
        pl.DataFrame({'sequence_id':ids,'clone_row':[0,0]}).write_csv(str(prefix)+'.read_ids.tsv',separator='\t')
        module.export_candidates(prefix,tmp_path/folder,representative_policy='native_root')
        chosen=pl.read_parquet(tmp_path/folder/'candidates.parquet').row(0,named=True)
        assert chosen['observed_read_id']=='alpha' and chosen['native_root_matching_observations']==2


def test_unknown_representative_policy_fails_before_output(tmp_path):
    prefix=native_fixture(tmp_path,assembled=False)
    with pytest.raises(ValueError,match='representative_policy'):
        module.export_candidates(prefix,tmp_path/'result',representative_policy='unsupported')
    assert not (tmp_path/'result').exists()


def test_native_root_short_observed_frame_retains_original_parent_search_entry(tmp_path):
    prefix=native_fixture(tmp_path,assembled=False)
    clones=pl.read_csv(str(prefix)+'.clones.tsv',separator='\t',infer_schema_length=0).to_dicts()
    clones[0].update(junction='TGTGCTTGG',junction_aa='CAW')
    pl.DataFrame(clones).write_csv(str(prefix)+'.clones.tsv',separator='\t')
    base=dict(rev_comp='F',locus='IGH',v_call='IGHV1-2',j_call='IGHJ4',c_call='IGHG1')
    root=dict(**base,sequence_id='short_root',sequence='TGTGCTTGG',junction='TGTGCTTGG',
        junction_aa='CAW',fwr3_start='1',fwr3_end='3',cdr3_start='4',cdr3_end='6',
        fwr4_start='7',fwr4_end='9',mmseqs2_score='100')
    wrong=dict(**base,sequence_id='long_wrong',sequence='GCT'*3+'TGTGTTTGG'+'GGT'*3,
        junction='TGTGTTTGG',junction_aa='CVW',fwr3_start='1',fwr3_end='12',
        cdr3_start='13',cdr3_end='15',fwr4_start='16',fwr4_end='27',mmseqs2_score='500')
    pl.DataFrame([root,wrong]).write_csv(str(prefix)+'.airr.tsv',separator='\t')
    pl.DataFrame({'sequence_id':['short_root','long_wrong'],'clone_row':[0,0]}).write_csv(str(prefix)+'.read_ids.tsv',separator='\t')
    old=module.export_candidates(prefix,tmp_path/'old')
    new=module.export_candidates(prefix,tmp_path/'new',representative_policy='native_root')
    assert old['counts']==new['counts']=={'catalog':3,'candidates':1,'search_entries':1}
    entry=pl.read_parquet(tmp_path/'new/search_entries.parquet').row(0,named=True)
    assert entry['seq']=='CAW' and entry['nt_seq']=='TGTGCTTGG'
    assert entry['observed_read_id']=='short_root' and entry['background_member'] is True
    # An incomplete long frame must NOT create search eligibility when the
    # original complete representative and corrected junction are both short.
    partial={**wrong,'sequence_id':'long_partial','sequence':'GCT'*3+'GCTGCTTGG'+'GGT'*3,
        'junction':'GCTGCTTGG','junction_aa':'AAW'}
    pl.DataFrame([root,partial]).write_csv(str(prefix)+'.airr.tsv',separator='\t')
    pl.DataFrame({'sequence_id':['short_root','long_partial'],'clone_row':[0,0]}).write_csv(str(prefix)+'.read_ids.tsv',separator='\t')
    old=module.export_candidates(prefix,tmp_path/'old_partial')
    new=module.export_candidates(prefix,tmp_path/'new_partial',representative_policy='native_root')
    assert old['counts']==new['counts']=={'catalog':3,'candidates':1,'search_entries':0}
