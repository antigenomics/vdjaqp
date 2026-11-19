"""Constructed source/native-PSM schemas; no cohort or native search execution."""
import gzip,json
from pathlib import Path
import polars as pl
import pytest
import vdjaqp.simulation_metrics as simulation_metrics
from scripts import export_all_locus_capture as report


def source_fixture(tmp_path,geometry='amplicon'):
    base=tmp_path/'observed';base.mkdir()
    seq='AAACDEFGHIKWAAA';junction='CDEFGHIKW'
    clones=pl.DataFrame({'locus':['IGH','IGH','IGK','IGL','IGH'],
        'junction_aa':[junction,junction,junction,'CAW','invalid'],
        'v_call':['IGHV1-2*01','IGHV1-2*01,IGHV3-23*01','IGKV1-5*01','IGLV1-40*01','IGHV1-2*01'],
        'j_call':['IGHJ4*01','IGHJ4*02','IGKJ1*01','IGLJ1*01','IGHJ4*01'],
        'c_call':['IGHG1*01','IGHG1*02',None,'IGLC1*01',None],
        'duplicate_count':[3,7,5,2,11]})
    clonefile=tmp_path/'smoke.clones.tsv';clones.write_csv(clonefile,separator='\t')
    rows=[];entries=[]
    for i,locus in enumerate(['IGH','IGH','IGK','IGL']):
        common=dict(source_observation_id=f'{geometry}_clone_{i}_read_0',parent_clone_id=f'{geometry}_clone_{i}',native_clone_row=i,locus=locus,native_junction_aa=junction if i<3 else 'CAW',duplicate_count=[3,7,5,2][i],background_member=True,coverage_kind='observed_variable_fragment' if i<3 else 'junction_only',sequence_origin='arda_observed_read' if i<3 else 'arda_corrected_junction')
        rows.append(common)
        if i<3:entries.append(dict(entry_id=f'entry_{i}',seq=seq,**common))
    frames={'catalog':pl.DataFrame(rows),'search_entries':pl.DataFrame(entries),'candidates':pl.DataFrame({'clone_id':['heavy0','heavy1'],'background_member':[True,True]})}
    for name,frame in frames.items():frame.write_parquet(base/f'{name}.parquet')
    receipt=dict(status='complete',geometry=geometry,source_format='arda_fastq_native',export_source_commit='53f1fd359303b2f76229058b325927c13556d655',counts={k:f.height for k,f in frames.items()},source_input_sha256={str(clonefile):report.digest(clonefile)},artifact_sha256={p.name:report.digest(p) for p in base.glob('*.parquet')},preprocess=dict(input_pairs=20,merged_pairs=12,unmerged_pairs=8,retained_output_records=28))
    path=base/'receipt.json';path.write_text(json.dumps(receipt))
    return path,clonefile


def acquisition_fixture(tmp_path,observed,geometry='amplicon'):
    paths=[]
    for block,seq in zip(report.BLOCKS,['AACDEF','CCCCCC','FGHIKW']):
        base=tmp_path/block;(base/'native').mkdir(parents=True)
        # Duplicate shared peptide assignment tests binary union and exact one-read acquisition.
        hits=''.join(f'<PeptideHit sequence="{seq}" charge="2" score="30"/>' for _ in range(2 if block==report.BLOCKS[2] else 1))
        xml=f'<IdXML><IdentificationRun><PeptideIdentification RT="1" MZ="500" score_type="hyperscore"><UserParam name="spectrum_reference" value="scan=1"/>{hits}</PeptideIdentification></IdentificationRun></IdXML>'
        native=base/'native/psms.idXML.gz';native.write_bytes(gzip.compress(xml.encode()))
        m=dict(status='complete',stage='acquisition',geometry=geometry,block=block,run_id=f'{geometry}:{block}',physical_H=1,physical_L=1 if block==report.BLOCKS[0] else 0,input_sha256={str(observed.resolve()):report.digest(observed)},artifact_sha256={'native/psms.idXML.gz':report.digest(native)})
        p=base/'manifest.json';p.write_text(json.dumps(m));paths.append(p)
    return paths


def fake_search_spy(monkeypatch):
    batches=[]
    def exact(targets,queries,threads):
        batches.append((targets.clone(),queries.clone(),threads));rows=[]
        for target,seq in targets.iter_rows():
            seq=seq.replace('I','L')
            for query,aa in queries.select('peptide_id','seq').iter_rows():
                aa=aa.replace('I','L');start=seq.find(aa)
                while start>=0:
                    rows.append((target,query,start,start+len(aa)));start=seq.find(aa,start+1)
        return pl.DataFrame(rows,schema={'chain_id':pl.String,'peptide_id':pl.String,'start':pl.Int32,'end':pl.Int32},orient='row')
    monkeypatch.setattr(simulation_metrics,'exact_peptide_placements',exact)
    return batches


def run_capture(tmp_path,monkeypatch):
    receipt,clones=source_fixture(tmp_path);paths=acquisition_fixture(tmp_path,receipt)
    _,records,inventory=report.load_observed(receipt,clones,'amplicon')
    calls,acquisitions,_=report.acquisition_calls(paths,'amplicon',receipt)
    batches=fake_search_spy(monkeypatch)
    tables,limits=report.capture(records,calls,cpus=2,max_placements=100,max_expanded_slots=10000)
    return tables,limits,batches,inventory,acquisitions


def test_duplicate_aliases_split_junction_class_unions_and_weights(tmp_path,monkeypatch):
    tables,limits,batches,inventory,acquisitions=run_capture(tmp_path,monkeypatch)
    assert len(batches)==limits['native_search_batches']==1
    assert batches[0][0].height==2 # one shared H/K segment plus short L junction target
    assert tables['records'].height==4 and tables['target_aliases'].height==4
    assert inventory.filter(pl.col('locus')=='IGH')['excluded_native_clone_records'].item()==1
    assert acquisitions.height==3 and acquisitions['psm_assignments'].sum()==4
    shared=acquisitions.filter(pl.col('acquisition')=='shared_chymotrypsin')
    assert shared['used_in_solution'].item() and shared['used_in_heavy_band'].item()
    c=tables['junction_class_coverage'].filter((pl.col('locus')=='IGH')&(pl.col('region')=='JUNCTION'))
    assert c.filter(pl.col('condition')=='solution')['coverage'].item()==1
    assert c.filter(pl.col('condition')=='heavy_band')['coverage'].item()==pytest.approx(6/9)
    assert c['member_records'].to_list()==[2,2] and c['native_support'].to_list()==[10,10]
    assert c.filter(pl.col('condition')=='solution')['complete_capture'].item()
    assert c.filter(pl.col('condition')=='heavy_band')['partial_only_capture'].item()
    short=tables['record_region_coverage'].filter(pl.col('locus')=='IGL')
    assert short.height==6 and short['covered_residues'].sum()==0 and not short['search_entry_present'].any()
    summary=tables['record_capture'].filter((pl.col('locus')=='IGH')&(pl.col('region')=='JUNCTION')&(pl.col('condition')=='solution'))
    assert summary['native_support_denominator'].item()==10 and summary['captured_native_support'].item()==10
    assert summary['weight_unit'].item()=='native_processed_fragment_support'
    ambiguity=tables['mapping_ambiguity'].filter(pl.col('compatible_native_records')>0)
    assert set(ambiguity['compatible_native_records'])=={3}
    assert set(ambiguity['compatible_exact_junction_classes'])=={2} # H and K stay distinct


def test_missing_constant_compound_calls_and_gene_denominators(tmp_path,monkeypatch):
    tables,*_=run_capture(tmp_path,monkeypatch)
    genes=tables['genes'];summary=tables['gene_capture']
    missing=genes.filter((pl.col('locus')=='IGK')&(pl.col('segment')=='C'))
    assert missing['call_status'].to_list()==['missing','missing']
    assert missing['native_support_denominator'].to_list()==[5,5]
    h=summary.filter((pl.col('locus')=='IGH')&(pl.col('segment')=='V')&(pl.col('condition')=='solution'))
    assert h['distinct_called_bins'].item()==2 and h['distinct_single_call_bins'].item()==1
    assert h['compound_or_unresolved_bins'].item()==1 and h['associated_called_bin_fraction'].item()==1
    j=summary.filter((pl.col('locus')=='IGH')&(pl.col('segment')=='J')&(pl.col('condition')=='solution'))
    assert j['distinct_called_bins'].item()==1 # allele suffixes share one unambiguous gene bin
    heavy_c=summary.filter((pl.col('locus')=='IGH')&(pl.col('segment')=='C')&(pl.col('condition')=='solution'))
    assert heavy_c['distinct_single_call_bins'].item()==1 # IGHG1 is a C call, not an IGHC-prefixed gene
    c=summary.filter((pl.col('locus')=='IGK')&(pl.col('segment')=='C'))
    assert c['distinct_called_bins'].to_list()==[0,0] and c['associated_called_bin_fraction'].null_count()==2


def test_il_equivalent_exact_junction_classes_stay_distinct(tmp_path,monkeypatch):
    receipt,clones=source_fixture(tmp_path)
    _,records,_=report.load_observed(receipt,clones,'amplicon')
    records=records.with_columns(pl.when(pl.col('native_clone_row')==1).then(pl.col('seq').str.replace_all('I','L')).otherwise(pl.col('seq')).alias('seq'),pl.when(pl.col('native_clone_row')==1).then(pl.col('native_junction_aa').str.replace_all('I','L')).otherwise(pl.col('native_junction_aa')).alias('native_junction_aa'))
    paths=acquisition_fixture(tmp_path,receipt);calls,_,_=report.acquisition_calls(paths,'amplicon',receipt);batches=fake_search_spy(monkeypatch)
    tables,_=report.capture(records,calls,cpus=1,max_placements=100,max_expanded_slots=10000)
    c=tables['junction_class_coverage'].filter((pl.col('locus')=='IGH')&(pl.col('condition')=='solution')&(pl.col('region')=='JUNCTION'))
    assert c.height==2 and c['native_junction_aa'].n_unique()==2 and c['complete_capture'].all()
    assert len(batches)==1


@pytest.mark.parametrize('bound',[1,3])
def test_native_hit_and_alias_expansion_limits(tmp_path,monkeypatch,bound):
    receipt,clones=source_fixture(tmp_path);_,records,_=report.load_observed(receipt,clones,'amplicon')
    paths=acquisition_fixture(tmp_path,receipt);calls,_,_=report.acquisition_calls(paths,'amplicon',receipt);fake_search_spy(monkeypatch)
    with pytest.raises(ValueError,match='(native hit|alias/condition)'):
        report.capture(records,calls,cpus=1,max_placements=bound,max_expanded_slots=10000)


def test_expansion_slot_limit(tmp_path,monkeypatch):
    receipt,clones=source_fixture(tmp_path);_,records,_=report.load_observed(receipt,clones,'amplicon')
    paths=acquisition_fixture(tmp_path,receipt);calls,_,_=report.acquisition_calls(paths,'amplicon',receipt);fake_search_spy(monkeypatch)
    with pytest.raises(ValueError,match='coverage expansion'):
        report.capture(records,calls,cpus=1,max_placements=100,max_expanded_slots=1)


def test_mixed_geometry_or_duplicate_acquisition_rejected(tmp_path):
    receipt,_=source_fixture(tmp_path);paths=acquisition_fixture(tmp_path,receipt)
    with pytest.raises(ValueError,match='geometry'):report.acquisition_calls(paths,'rnaseq',receipt)
    with pytest.raises(ValueError,match='distinct'):report.acquisition_calls([paths[0],paths[0],paths[2]],'amplicon',receipt)


def test_native_clone_source_hash_mismatch_rejected(tmp_path):
    receipt,clones=source_fixture(tmp_path);clones.write_text(clones.read_text().replace('IGHG1*01','IGHG2*01'))
    with pytest.raises(ValueError,match='source hash'):report.load_observed(receipt,clones,'amplicon')


def test_full_cli_units_hashes_fresh_output_and_one_batch(tmp_path,monkeypatch):
    receipt,clones=source_fixture(tmp_path);paths=acquisition_fixture(tmp_path,receipt);batches=fake_search_spy(monkeypatch)
    output=tmp_path/'report';args=['--observed-receipt',str(receipt),'--native-clones',str(clones),'--geometry','amplicon','--cpus','1','--output',str(output)]
    for path in paths:args+=['--acquisition-manifest',str(path)]
    report.main(args);m=json.loads((output/'manifest.json').read_text())
    assert len(batches)==1 and m['limits']['native_search_batches']==1
    totals=pl.read_parquet(output/'input_totals.parquet')
    assert totals['count'].to_list()==[20,28,17]
    assert totals['unit'].to_list()==['raw_read_pairs','processed_fragment_records','native_processed_fragment_support']
    assert all(report.digest(output/name)==h for name,h in m['artifact_sha256'].items())
    assert all(report.digest(name)==h for name,h in m['source_input_sha256'].items())
    with pytest.raises(FileExistsError):report.main(args)


@pytest.mark.parametrize('change',['missing_long_entry','inserted_role','wrong_entry_provenance'])
def test_invalid_observed_contract_rejected(tmp_path,change):
    receipt,clones=source_fixture(tmp_path);m=json.loads(receipt.read_text());base=receipt.parent
    name='search_entries' if change!='inserted_role' else 'catalog'
    f=pl.read_parquet(base/f'{name}.parquet')
    if change=='missing_long_entry':f=f.slice(1)
    elif change=='inserted_role':f=f.with_columns(background_member=pl.lit(False))
    else:f=f.with_columns(parent_clone_id=pl.lit('wrong_parent'))
    f.write_parquet(base/f'{name}.parquet');m['artifact_sha256'][name+'.parquet']=report.digest(base/f'{name}.parquet');m['counts'][name]=f.height;receipt.write_text(json.dumps(m))
    with pytest.raises(ValueError,match='(missing search|background-role|provenance mismatch)'):
        report.load_observed(receipt,clones,'amplicon')


def test_shared_acquisition_read_once_and_pool_aliases(tmp_path,monkeypatch):
    receipt,_=source_fixture(tmp_path);paths=acquisition_fixture(tmp_path,receipt);seen=[]
    import vdjaqp.readout
    reader=vdjaqp.readout.read_idxml
    def read(path,*,run_id):seen.append(run_id);return reader(path,run_id=run_id)
    monkeypatch.setattr(vdjaqp.readout,'read_idxml',read)
    calls,inventory,_=report.acquisition_calls(paths,'amplicon',receipt)
    assert len(seen)==len(set(seen))==3 and sum(':shared_chymotrypsin' in x for x in seen)==1
    assert calls.filter(pl.col('condition')=='solution').height==3
    assert calls.filter(pl.col('condition')=='heavy_band').height==3
    assert inventory['psm_assignments'].sum()==4


def test_zero_hits_preserve_all_native_denominators(tmp_path,monkeypatch):
    receipt,clones=source_fixture(tmp_path);_,records,_=report.load_observed(receipt,clones,'amplicon')
    paths=acquisition_fixture(tmp_path,receipt);calls,_,_=report.acquisition_calls(paths,'amplicon',receipt)
    calls=calls.with_columns(seq=pl.lit('PPPPPP'));fake_search_spy(monkeypatch)
    tables,limits=report.capture(records,calls,cpus=1,max_placements=100,max_expanded_slots=10000)
    assert limits['native_hits']==0 and tables['record_region_coverage'].height==24
    assert tables['record_region_coverage']['covered_residues'].sum()==0
    assert tables['mapping_ambiguity']['compatible_native_records'].sum()==0
    assert tables['record_capture'].filter((pl.col('condition')=='solution')&(pl.col('region')=='JUNCTION'))['native_support_denominator'].sum()==17


def test_short_only_queries_do_not_create_matching_or_denominator_loss(tmp_path,monkeypatch):
    receipt,clones=source_fixture(tmp_path);_,records,_=report.load_observed(receipt,clones,'amplicon')
    paths=acquisition_fixture(tmp_path,receipt);calls,_,_=report.acquisition_calls(paths,'amplicon',receipt)
    calls=calls.with_columns(seq=pl.lit('CAW'));fake_search_spy(monkeypatch)
    tables,limits=report.capture(records,calls,cpus=1,max_placements=100,max_expanded_slots=10000)
    assert tables['queries'].height==0 and tables['excluded_short_queries'].height==3
    assert limits['native_search_batches']==0 and tables['record_region_coverage'].height==24
    assert tables['record_region_coverage']['covered_residues'].sum()==0


def test_slash_separated_compound_gene_call_not_single_call(tmp_path,monkeypatch):
    receipt,clones=source_fixture(tmp_path);_,records,_=report.load_observed(receipt,clones,'amplicon')
    records=records.with_columns(pl.when(pl.col('native_clone_row')==1).then(pl.lit('IGHV1-2/IGHV3-23')).otherwise(pl.col('v_call')).alias('v_call'))
    paths=acquisition_fixture(tmp_path,receipt);calls,_,_=report.acquisition_calls(paths,'amplicon',receipt);fake_search_spy(monkeypatch)
    tables,_=report.capture(records,calls,cpus=1,max_placements=100,max_expanded_slots=10000)
    compound=tables['genes'].filter(pl.col('call_bin')=='IGHV1-2/IGHV3-23')
    assert compound.height==2 and set(compound['call_status'])=={'compound_or_unresolved'}


def test_all_short_catalog_rows_survive_empty_search_entry_table(tmp_path):
    receipt,clones=source_fixture(tmp_path);base=receipt.parent;m=json.loads(receipt.read_text())
    native=pl.read_csv(clones,separator='\t').with_columns(junction_aa=pl.lit('CAW'))
    native.write_csv(clones,separator='\t');m['source_input_sha256'][str(clones)]=report.digest(clones)
    catalog=pl.read_parquet(base/'catalog.parquet').with_columns(native_junction_aa=pl.lit('CAW'),coverage_kind=pl.lit('junction_only'))
    entries=pl.DataFrame(schema={'entry_id':pl.String,'seq':pl.String,'source_observation_id':pl.String,'background_member':pl.Boolean})
    for name,f in [('catalog',catalog),('search_entries',entries)]:
        f.write_parquet(base/f'{name}.parquet');m['artifact_sha256'][name+'.parquet']=report.digest(base/f'{name}.parquet');m['counts'][name]=f.height
    receipt.write_text(json.dumps(m));_,records,_=report.load_observed(receipt,clones,'amplicon')
    assert records.height==4 and not records['search_entry_present'].any() and records['seq'].to_list()==['CAW']*4


def test_direct_canonical_producer_without_private_commit_is_accepted(tmp_path,monkeypatch):
    receipt,clones=source_fixture(tmp_path);m=json.loads(receipt.read_text());m.pop('export_source_commit');receipt.write_text(json.dumps(m))
    _,records,_=report.load_observed(receipt,clones,'amplicon');assert records.height==4
    paths=acquisition_fixture(tmp_path,receipt);fake_search_spy(monkeypatch)
    output=tmp_path/'direct_report';args=['--observed-receipt',str(receipt),'--native-clones',str(clones),'--geometry','amplicon','--cpus','1','--output',str(output)]
    for path in paths:args+=['--acquisition-manifest',str(path)]
    report.main(args);manifest=json.loads((output/'manifest.json').read_text())
    assert manifest['status']=='complete' and manifest['observed_producer_provenance']=={}


@pytest.mark.parametrize('change',['source_format','catalog_field','entry_field'])
def test_wrong_format_or_missing_required_fields_rejected(tmp_path,change):
    receipt,clones=source_fixture(tmp_path);m=json.loads(receipt.read_text())
    if change=='source_format':m['source_format']='not_native_arda'
    else:
        name,column=('catalog','native_clone_row') if change=='catalog_field' else ('search_entries','seq')
        path=receipt.parent/f'{name}.parquet';pl.read_parquet(path).drop(column).write_parquet(path)
        m['artifact_sha256'][path.name]=report.digest(path)
    receipt.write_text(json.dumps(m))
    with pytest.raises(ValueError,match='(matching native ARDA|required.*schema)'):
        report.load_observed(receipt,clones,'amplicon')



def test_cli_module_load_precedes_numerical_imports():
    import subprocess,sys
    subprocess.run([sys.executable,'-c',"import runpy,sys;runpy.run_path('scripts/export_all_locus_capture.py',run_name='import_check');assert 'polars' not in sys.modules and 'numpy' not in sys.modules"],check=True,timeout=10)


def test_native_c_class_remains_call_bin_without_gene_claim(tmp_path,monkeypatch):
    receipt,clones=source_fixture(tmp_path);_,records,_=report.load_observed(receipt,clones,'amplicon')
    records=records.with_columns(pl.when(pl.col('locus')=='IGH').then(pl.lit('IGHG')).otherwise(pl.col('c_call')).alias('c_call'))
    paths=acquisition_fixture(tmp_path,receipt);calls,_,_=report.acquisition_calls(paths,'amplicon',receipt);fake_search_spy(monkeypatch)
    tables,_=report.capture(records,calls,cpus=1,max_placements=100,max_expanded_slots=10000)
    c=tables['genes'].filter((pl.col('locus')=='IGH')&(pl.col('segment')=='C'))
    assert set(c['call_bin'])=={'IGHG'} and set(c['call_status'])=={'single_call'}
    assert set(c['call_namespace'])=={'native_C_call_not_gene_or_subclass_verified'}
    summary=tables['gene_capture']
    assert 'distinct_single_call_bins' in summary.columns and 'distinct_single_gene_bins' not in summary.columns
