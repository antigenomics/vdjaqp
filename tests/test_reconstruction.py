"""Known-reference insertion is explicit and cannot fabricate RNA support."""
import polars as pl
import pytest
from vdjaqp.reconstruction import append_known_references


def inputs():
    nt='GCT'*6+'TGTGCTTTT'+'GCT'*6
    candidate=pl.DataFrame([dict(clone_id='s_c',seq='AAAAAACAF'+'AAAAAA',nt_seq=nt,
        cdr3_seq='CAF',cdr3_start=6,cdr3_end=9,region_ids=[4]*7+[5]+[7]*7,
        v_name='IGHV1-2*01',j_name='IGHJ4*01',c_name='IGHG1',background_member=True,
        sequence_origin='arda_observed_read',coverage_kind='observed_variable_fragment',
        native_junction_nt='TGTGCTTTT',native_junction_aa='CAF',native_clone_row=0,
        parent_clone_id='native0',observed_read_id='read0')])
    entry=candidate.select(pl.col('clone_id').alias('entry_id'),'seq','background_member',
        'parent_clone_id','sequence_origin').with_columns(pl.lit('IGH').alias('locus'),
        pl.lit('read0').alias('source_observation_id'))
    observed={'candidates':candidate,'search_entries':entry}
    reference=observed['candidates'].select(
        'clone_id','seq','nt_seq','cdr3_seq','cdr3_start','cdr3_end',
        'region_ids','v_name','j_name','c_name').with_columns(
        pl.lit('known0').alias('clone_id'),pl.lit('IGHG1').alias('c_name'))
    return observed['candidates'],observed['search_entries'],reference


@pytest.mark.parametrize('count',[1,12])
def test_appends_without_relabelling_observed_support(count):
    candidates,entries,reference=inputs()
    reference=pl.concat([reference.with_columns(pl.lit(f'known{i}').alias('clone_id')) for i in range(count)])
    result=append_known_references(candidates,entries,reference,
        expected_count=count,source_format='synthetic known-panel annotation')
    c,s=result['candidates'],result['search_entries']
    assert c.height==1+count and s.height==entries.height+count
    assert c.filter(pl.col('clone_id')=='s_c').select(candidates.columns).equals(candidates)
    inserted=c.filter(pl.col('clone_id')=='known0').row(0,named=True)
    assert inserted['candidate_construction']=='reference_panel_insertion'
    assert inserted['coverage_kind']=='supplied_reference_sequence_not_rna_support'
    assert inserted['observed_read_id'] is None
    assert c['background_member'].to_list()==[True]+[False]*count
    assert s['background_member'].to_list()==[True]+[False]*count
    assert result['reference_panel']['igseq_source_format'][0]=='synthetic known-panel annotation'
    assert s.filter(pl.col('entry_id')=='known0')['sequence_origin'][0]=='supplied_known_reference'


def test_rejects_collision_wrong_count_and_invalid_sequence():
    c,s,p=inputs()
    with pytest.raises(ValueError,match='count'):
        append_known_references(c,s,p,expected_count=12,source_format='test')
    for key,value in [('clone_id','s_c'),('nt_seq','AAA'),('v_name','IGKV1-1'),('v_name',None)]:
        with pytest.raises(ValueError):
            append_known_references(c,s,p.with_columns(pl.lit(value).alias(key)),
                expected_count=1,source_format='test')
    first=append_known_references(c,s,p,expected_count=1,source_format='test')
    with pytest.raises(ValueError):
        append_known_references(first['candidates'],first['search_entries'],p,
            expected_count=1,source_format='test')


def test_cli_records_counts_hashes_and_never_overwrites(tmp_path):
    import hashlib,json,subprocess,sys
    c,s,p=inputs()
    for name,frame in [('c',c),('s',s),('p',p)]:frame.write_parquet(tmp_path/(name+'.parquet'))
    out=tmp_path/'new'
    command=[sys.executable,'-m','vdjaqp.reconstruction','--candidates',str(tmp_path/'c.parquet'),
        '--search-entries',str(tmp_path/'s.parquet'),'--reference-panel',str(tmp_path/'p.parquet'),
        '--expected-count','1','--source-format','synthetic fixture','--output',str(out)]
    subprocess.run(command,check=True,capture_output=True,timeout=30)
    m=json.loads((out/'insertion_manifest.json').read_text())
    assert m['input_rows']=={'candidates':1,'search_entries':1,'reference_panel':1}
    assert m['output_rows']['candidates']==2
    assert m['input_sha256']['reference_panel']==hashlib.sha256((tmp_path/'p.parquet').read_bytes()).hexdigest()
    original=(out/'candidates.parquet').read_bytes()
    assert subprocess.run(command,capture_output=True,timeout=30).returncode!=0
    assert (out/'candidates.parquet').read_bytes()==original


def test_insertion_preserves_explicit_background_role_and_rejects_missing():
    c,s,p=inputs()
    for value in [True,False]:
        tagged=c.with_columns(background_member=pl.lit(value))
        out=append_known_references(tagged,s,p,expected_count=1,source_format='fixture')['candidates']
        assert out['background_member'].to_list()==[value,False]
        assert out.head(1).select(tagged.columns).equals(tagged)
    for value in [None,1,'true']:
        with pytest.raises(ValueError,match='background_member'):
            append_known_references(c.with_columns(background_member=pl.lit(value)),s,p,
                                    expected_count=1,source_format='fixture')


@pytest.mark.parametrize('origin',['arda_observed_read','arda_native_assembled_read','arda_corrected_junction'])
def test_native_origin_and_coding_provenance(origin):
    c,s,p=inputs()
    if origin=='arda_corrected_junction':
        c=c.with_columns(seq=pl.lit('CAF'),nt_seq=pl.lit('TGTGCTTTT'),
            cdr3_start=pl.lit(0),cdr3_end=pl.lit(3),region_ids=pl.lit([4,5,7]))
    c=c.with_columns(sequence_origin=pl.lit(origin),coverage_kind=pl.lit(
        'junction_only' if origin=='arda_corrected_junction' else 'observed_variable_fragment'))
    out=append_known_references(c,s,p,expected_count=1,source_format='fixture')
    assert out['candidates'].head(1).select(c.columns).equals(c)
    assert out['search_entries'].head(1).select(s.columns).equals(s)
    for broken in [c.drop('native_junction_nt'),c.with_columns(nt_seq=pl.lit('AAA')),
                   c.with_columns(sequence_origin=pl.lit('native_preimputation_contig_retained_primary_support')),
                   c.with_columns(native_junction_aa=pl.lit('CCF'))]:
        with pytest.raises(ValueError,match='native ARDA'):
            append_known_references(broken,s,p,expected_count=1,source_format='fixture')
