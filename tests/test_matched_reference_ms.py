"""Constructed bridge tests; no native annotation, spectra or scientific cohort."""
import polars as pl
import pytest
from scripts.run_matched_reference_ms import aa_regions, class_benchmark, acquisition_scope, build_database, retire_native



def arda_inputs():
    """Native-origin coding fragments, no assembly/mapping fixture dependency."""
    candidate=pl.DataFrame({'clone_id':['observed0'],'seq':['ACAFWG'],'nt_seq':['GCTTGTGCTTTTTGGGGT'],'cdr3_seq':['CAFW'],'cdr3_start':[1],'cdr3_end':[5],'region_ids':[[0,4,6,6,7,7]],'v_name':['IGHV1-2*01'],'j_name':['IGHJ4*02'],'c_name':['IGHM'],'sequence_origin':['arda_observed_read'],'coverage_kind':['observed_variable_fragment'],'native_clone_row':[0],'parent_clone_id':['native0'],'source_observation_id':['native_read0'],'duplicate_count':[3],'native_junction_nt':['TGTGCTTTTTGG'],'native_junction_aa':['CAFW'],'background_member':[True]})
    entries=candidate.select(pl.lit('entry0').alias('entry_id'),'seq','source_observation_id','parent_clone_id','sequence_origin','coverage_kind','native_clone_row','native_junction_nt','native_junction_aa').with_columns(pl.lit('IGH').alias('locus'),pl.lit(True).alias('background_member'))
    panel=candidate.select('clone_id','seq','nt_seq','cdr3_seq','cdr3_start','cdr3_end','region_ids','v_name','j_name','c_name').with_columns(pl.lit('known0').alias('clone_id'),pl.lit('IGHG1').alias('c_name'))
    return candidate,entries,panel


def test_preparation_inventory_is_three_unique_acquisitions():
    assert acquisition_scope('solution_trypsin')==('Trypsin',False)
    assert acquisition_scope('heavy_band_trypsin')==('Trypsin',True)
    assert acquisition_scope('shared_chymotrypsin')==('Chymotrypsin',True)
    with pytest.raises(ValueError):acquisition_scope('new_condition')


def test_aa_units_and_no_invented_nt_d_or_n():
    m=pl.DataFrame({'chain_id':['h'],'locus':['IGH'],'variable_aa':['ACDEFGW'],'protein_aa':['ACDEFGWAAA'],'variable_length':[7],'known_junction':['CDEFGW']})
    a=[{'sequence_id':'h','sequence':'ACDEFGW','fwr1_start':1,'fwr1_end':2,'fwr1_aa':'AC','v_call':'V','j_call':'J'}]
    r=aa_regions(m,a)
    assert r.filter(pl.col('region')=='FR1').select('start','end').row(0)==(0,2)
    assert r.filter(pl.col('region')=='JUNCTION').select('start','end').row(0)==(1,7)
    assert r.filter(pl.col('region').is_in(['D','NP1','NP2']))['start'].null_count()==3
    assert r.filter(pl.col('region')=='C').select('start','end').row(0)==(7,10)


def test_class_universe_aliases_ambiguity_missing_sentinel_and_unresolved_calls():
    universe=pl.DataFrame({'variable_mass_class_id':['pos','pos','neg','missing'],'physical_positive':[True,False,False,False]})
    pairs=pl.DataFrame({'clone_id':['q','alias','negative','ambiguous','ambiguous'],'variable_mass_class_id':['pos','pos','neg','pos','neg']})
    scores=pl.DataFrame({'clone_id':['q','alias','negative','ambiguous','partial'],'score':[-5.,-6.,-10.,50.,70.],'confirmed_nats':[21.,22.,25.,30.,30.]})
    classes,mapping,metrics=class_benchmark(universe,pairs,scores)
    assert classes.height==3
    assert classes.filter(pl.col('clone_id')=='missing')['score'].item()==-11
    assert classes.filter(pl.col('clone_id')=='pos')['score'].item()==-5
    assert classes.filter(pl.col('clone_id')=='pos')['confirmed_nats'].item()==22
    assert mapping.filter(pl.col('clone_id')=='ambiguous')['full_class_count'].item()==2
    assert mapping.filter(pl.col('clone_id')=='ambiguous')['variable_mass_class_id'].item() is None
    assert metrics['all_shortlist_candidate_records']['unresolved_called_records']==2
    assert metrics['all_shortlist_candidate_records']['definitive_precision_all_calls']==2/5
    assert metrics['end_to_end_full_variable_IL_classes']['n_positives']==1
    assert metrics['database_observable_full_variable_IL_classes']['n_candidates']==2


def test_database_never_accepts_full_physical_model(tmp_path):
    entries=pl.DataFrame({'entry_id':['observed','query'],'seq':['ACDEFG','HIKLMN']})
    # Explicit known-source auxiliaries/constants, never an argument for physical proteins.
    aux={'known_light':'NPQRST'};constants={'public_C':'VWYACD'}
    p=tmp_path/'db.fasta';build_database(entries,aux,constants,p)
    assert 'PHYSICAL' not in p.read_text() and p.read_text().count('>')==4
    with pytest.raises(ValueError):build_database(entries,{'query':'NPQRST'},constants,tmp_path/'collision')


def test_retirement_only_after_complete_raw_export_validation(tmp_path):
    native=tmp_path/'native';native.mkdir()
    for name in ['spectra.mzML','truth.featureXML','psms.idXML','peptides.tsv','ms2_truth.json']:(native/name).write_text(name)
    with pytest.raises(ValueError):retire_native(native,validated=False)
    ledger=retire_native(native,validated=True)
    assert len(ledger)==2 and not (native/'spectra.mzML').exists()
    assert (native/'psms.idXML').exists() and (native/'peptides.tsv').exists()


@pytest.mark.parametrize('block,expected_H,expected_L',[('solution_trypsin',12,12),('heavy_band_trypsin',12,0),('shared_chymotrypsin',12,0)])
def test_first_acquisition_path_only_uses_declared_sources(tmp_path,monkeypatch,block,expected_H,expected_L):
    import json,gzip
    from types import SimpleNamespace
    from scripts import run_matched_reference_ms as bridge
    from test_reconstruction import inputs
    import vdjaqp.readout
    import vdjaqp.simulation_metrics
    background,entries,one=inputs()
    background=background.with_columns(pl.lit(True).alias('background_member'))
    entries=entries.with_columns(pl.lit(True).alias('background_member'))
    panel=pl.concat([one.with_columns(pl.lit(f'known{i}').alias('clone_id')) for i in range(12)])
    packet=tmp_path/'packet';packet.mkdir();manifest=packet/'packet_manifest.json';manifest.write_text(json.dumps({'status':'complete','artifact_sha256':{}}))
    (packet/'supplied_light_auxiliary_search.fasta').write_text('>supplied_light\nACDEFG\n')
    (packet/'public_constant_search.fasta').write_text('>public_C\nHIKLMN\n')
    physical=pl.DataFrame([dict(chain_id=f'ref{i}_{side}',locus='IGH' if side=='H' else 'IGK',variable_aa='ACDEFG',protein_aa='ACDEFGAAA',variable_length=6) for i in range(12) for side in ['H','L']])
    obs=tmp_path/'observed';obs.mkdir();background.write_parquet(obs/'candidates.parquet');entries.write_parquet(obs/'search_entries.parquet')
    receipt=tmp_path/'receipt.json';receipt.write_text('{}')
    regions=tmp_path/'regions.json';regions.write_text(json.dumps({'status':'complete','artifact_sha256':{},'packet_manifest_sha256':bridge.digest(manifest)}))
    simulator=tmp_path/'simulator';simulator.write_text('constructed binary pin, never executed')
    monkeypatch.setattr(bridge,'packet_frames',lambda path:({},panel,physical))
    monkeypatch.setattr(bridge,'verify_observations',lambda *a:({}, {},obs))
    commands=[]
    def native(command,**kwargs):
        commands.append(command);dest=__import__('pathlib').Path(command[4]);dest.mkdir()
        assert '--keep-spectra' in command and command[command.index('--ms2-per-bin')+1]=='5'
        protein_input=__import__('pathlib').Path(command[2]).read_text();database=__import__('pathlib').Path(command[3]).read_text()
        assert protein_input.count('>')==expected_H+expected_L
        assert '[# intensity=10000 #]' in protein_input
        assert 'ACDEFGAAA' not in database # physical completions do not enter search DB
        for name in ['spectra.mzML','truth.featureXML','truth.idXML','peptides.tsv']:(dest/name).write_text(name)
        (dest/'psms.idXML').write_text('<IdXML><IdentificationRun><PeptideIdentification RT="10" MZ="500" score_type="hyperscore"><UserParam name="spectrum_reference" value="scan=9"/><PeptideHit sequence="ACDEFG" charge="2" score="33"/></PeptideIdentification></IdentificationRun></IdXML>')
        (dest/'ms2_truth.json').write_text(json.dumps([{'native_spectrum_reference':'scan=9','contributors':[{'peptide_aa':'ACDEFG','protein_accessions':['ref0_H']}]}]))
    monkeypatch.setattr(bridge.subprocess,'run',native)
    output=tmp_path/'out';args=SimpleNamespace(packet_manifest=manifest,reference_regions=regions,observed_receipt=receipt,geometry='amplicon',block=block,output=output,cpus=1,mssimulator=simulator,stage_timeout=60)
    bridge.acquire(args)
    assert len(commands)==1
    m=json.loads((output/'manifest.json').read_text());assert m['physical_H']==expected_H and m['physical_L']==expected_L
    c=pl.read_parquet(output/'candidates.parquet');assert c.height==background.height+12
    assert c.filter(pl.col('clone_id').str.starts_with('known'))['background_member'].to_list()==[False]*12
    assert not (output/'native/spectra.mzML').exists() and not (output/'native/truth.featureXML').exists()
    with gzip.open(output/'native/psms.idXML.gz','rb') as f:parsed=vdjaqp.readout.read_idxml(f,run_id=m['run_id'])
    assert parsed.height==1
    assert pl.read_parquet(output/'psm_accuracy.parquet').filter(pl.col('quantity')=='correct_psm_rate')['fraction'].item()==1


def test_one_geometry_uses_three_acquisitions_two_pools_one_exact_batch(tmp_path,monkeypatch):
    import gzip,json
    from types import SimpleNamespace
    from scripts import run_matched_reference_ms as bridge
    import vdjaqp.pipeline,vdjaqp.simulation_metrics,vdjaqp.airr_capture
    packet=tmp_path/'packet.json';packet.write_text('{}');outer=tmp_path/'outer.json';outer.write_text('{}')
    refdir=tmp_path/'region';refdir.mkdir();manifest=refdir/'manifest.json'
    metadata=pl.DataFrame({'chain_id':['ref_H','ref_L'],'locus':['IGH','IGK'],'variable_aa':['ACDEFG','HIKLMN'],'protein_aa':['ACDEFGAAA','HIKLMNAAA'],'variable_length':[6,6],'known_junction':[None,None],'constant_gene':['IGHG1','IGKC'],'constant_nt_gene':['IGHG1','IGKC'],'v_call':['IGHV1*01','IGKV1*01'],'j_call':['IGHJ1*01','IGKJ1*01'],'canonical_v_gene':[None,None],'canonical_v_match_status':['native_AA_inferred_not_canonical_mapped']*2,'canonical_v_gene_count':[None,None]})
    annotations=[dict(sequence_id=row['chain_id'],sequence=row['variable_aa'],fwr1_start=1,fwr1_end=6,fwr1_aa=row['variable_aa'],v_sequence_start=1,v_sequence_end=6,v_call=row['v_call'],j_call=row['j_call'])for row in metadata.iter_rows(named=True)]
    regions=aa_regions(metadata,annotations);metadata.write_parquet(refdir/'physical_metadata.parquet');regions.write_parquet(refdir/'physical_regions.parquet')
    manifest.write_text(json.dumps({'status':'complete','packet_manifest_sha256':bridge.digest(packet),'artifact_sha256':{p.name:bridge.digest(p)for p in refdir.glob('*.parquet')}}))
    candidates=pl.DataFrame({'clone_id':['known','background'],'seq':['ACDEFGAAA','CCCCCC'],'c_name':['IGHG1','IGHM'],'v_name':['IGHV1*01']*2,'j_name':['IGHJ1*01']*2,'background_member':[False,True]})
    monkeypatch.setattr(bridge,'inference_inputs',lambda *a:({'candidates':candidates,'reference_panel':pl.DataFrame({'sequence_id':['known']})},metadata))
    monkeypatch.setattr(bridge,'verify_observations',lambda *a:({}, {},tmp_path))
    bg=pl.DataFrame({'chain_id':['bg_H'],'locus':['IGH'],'variable_aa':['CCCCCC'],'protein_aa':['CCCCCCAAA'],'variable_length':[6]})
    monkeypatch.setattr(bridge,'raw_background',lambda *a:bg)
    truth=tmp_path/'truth.gz';truth.write_bytes(b'not read: helper mocked');validation=tmp_path/'validation.json';validation.write_text('{}')
    manifests=[]
    for block in bridge.BLOCKS:
        d=tmp_path/block;(d/'native').mkdir(parents=True);peptide='HIKLMN' if block==bridge.BLOCKS[0] else 'ACDEFG';protein='ref_L' if block==bridge.BLOCKS[0] else 'ref_H'
        xml=f'<IdXML><IdentificationRun><PeptideIdentification RT="10" MZ="500" score_type="hyperscore"><UserParam name="spectrum_reference" value="scan=9"/><PeptideHit sequence="{peptide}" charge="2" score="33"/></PeptideIdentification></IdentificationRun></IdXML>'
        (d/'native/psms.idXML.gz').write_bytes(gzip.compress(xml.encode()));(d/'native/ms2_truth.json.gz').write_bytes(gzip.compress(json.dumps([{'native_spectrum_reference':'scan=9','contributors':[{'peptide_aa':peptide,'protein_accessions':[protein]}]}]).encode()))
        m=dict(status='complete',stage='acquisition',geometry='amplicon',block=block,physical_H=1,physical_L=1 if block==bridge.BLOCKS[0] else 0,run_id='amplicon:'+block,input_sha256={str(p.resolve()):bridge.digest(p)for p in [packet,manifest,outer]},artifact_sha256={str(p.relative_to(d)):bridge.digest(p)for p in (d/'native').iterdir()})
        f=d/'manifest.json';f.write_text(json.dumps(m));manifests.append(f)
    evidence_calls=[]
    class FakeEvidence:
        @staticmethod
        def from_frames(clones,calls,error_model):
            assert set(clones['clone_id'])=={'known','background'} and 'protein_aa' not in clones.columns
            evidence_calls.append(calls['peptide_id'].to_list());return FakeEvidence()
        def score(self):return pl.DataFrame({'clone_id':['known','background'],'score':[30.,25.],'confirmed_nats':[30.,25.]})
    monkeypatch.setattr(vdjaqp.pipeline,'Evidence',FakeEvidence)
    monkeypatch.setattr(vdjaqp.airr_capture,'capture_tables',lambda *a:{})
    exact_calls=[]
    def exact(chains,queries,threads):
        exact_calls.append(queries.height);out=[]
        for chain,protein in chains.iter_rows():
            for seq,q in queries.select('seq','peptide_id').iter_rows():
                start=protein.find(seq)
                if start>=0:out.append((chain,q,start,start+len(seq)))
        return pl.DataFrame(out,schema={'chain_id':pl.String,'peptide_id':pl.String,'start':pl.Int32,'end':pl.Int32},orient='row')
    monkeypatch.setattr(vdjaqp.simulation_metrics,'exact_peptide_placements',exact)
    output=tmp_path/'evaluation'
    args=SimpleNamespace(packet_manifest=packet,reference_regions=manifest,observed_receipt=outer,geometry='amplicon',background_truth_gzip=truth,background_validation_receipt=validation,acquisition_manifests=manifests,cpus=1,max_placements=1000,max_expanded_slots=10000,output=output)
    bridge.evaluate(args)
    assert len(exact_calls)==1 and len(evidence_calls)==2
    assert set(evidence_calls[0])&set(evidence_calls[1])=={'amplicon:shared_chymotrypsin:psm:0'}
    coverage=pl.read_parquet(output/'heavy_band_physical_region_coverage.parquet').filter(pl.col('locus')=='IGK')
    assert coverage['coverage'].null_count()==coverage.height and set(coverage['status'])=={'not_assayed'}
    metrics=json.loads((output/'metrics.json').read_text());assert set(metrics)=={'solution','heavy_band'}
    assert metrics['solution']['end_to_end_full_variable_IL_classes']['n_positives']==1


def test_separate_source_region_phase_is_one_aa_batch(tmp_path,monkeypatch):
    import json
    from types import SimpleNamespace
    from scripts import run_matched_reference_ms as bridge
    import arda.annotate.contig
    metadata=pl.DataFrame({'chain_id':['h'],'locus':['IGH'],'variable_aa':['ACDEFGW'],'protein_aa':['ACDEFGWAAA'],'variable_length':[7],'known_junction':['CDEFGW'],'constant_gene':['IGHG1']})
    packet=tmp_path/'packet.json';packet.write_text('{}');monkeypatch.setattr(bridge,'packet_frames',lambda p:({},None,metadata))
    calls=[]
    def annotate(records,**kwargs):
        calls.append((records,kwargs));assert kwargs['seqtype']=='aa' and kwargs['map_d'] is False
        return [dict(sequence_id='h',sequence='ACDEFGW',v_call='IGHV1*01',j_call='IGHJ1*01',fwr1_start=1,fwr1_end=2,fwr1_aa='AC')]
    monkeypatch.setattr(arda.annotate.contig,'reannotate_contigs',annotate)
    output=tmp_path/'regions';bridge.qualify_regions(SimpleNamespace(packet_manifest=packet,output=output,cpus=1))
    assert len(calls)==1 and calls[0][0]==[('h','ACDEFGW')]
    assert json.loads((output/'manifest.json').read_text())['annotation_batches']==1
    assert pl.read_parquet(output/'physical_metadata.parquet')['canonical_v_match_status'].item()=='native_AA_inferred_not_canonical_mapped'


def test_cli_load_has_no_numeric_import_or_native_execution():
    import subprocess,sys
    subprocess.run([sys.executable,'-c',"import runpy,sys;runpy.run_path('scripts/run_matched_reference_ms.py',run_name='load_check');assert 'numpy' not in sys.modules and 'polars' not in sys.modules"],check=True,timeout=10)



def test_experimental_readouts_reuse_query_masks_hits_and_shared_acquisition(monkeypatch):
    from scripts.run_matched_reference_ms import experimental_readouts,BLOCKS
    from vdjaqp.airr_ms_metrics import supplied_reference_regions
    import vdjaqp.pipeline,vdjaqp.simulation_metrics
    def forbidden(*args,**kwargs):raise AssertionError('reporting must not search or score')
    monkeypatch.setattr(vdjaqp.pipeline,'Evidence',forbidden)
    monkeypatch.setattr(vdjaqp.simulation_metrics,'exact_peptide_placements',forbidden)
    variable='ALGTSCAAWFG'
    panel=pl.DataFrame([dict(clone_id=f'q{i}',seq=variable+'AAA',region_ids=[0,1,2,3,4,5,6,6,5,7,7,8,8,8],cdr3_seq='CAAW',cdr3_start=5,cdr3_end=9)for i in range(1,13)])
    physical=pl.DataFrame([dict(chain_id=f'h{i}',locus='IGH',query_clone_id=f'q{i}',reference=f'ref{i:02d}',variable_length=11,variable_aa=variable,protein_aa=variable+'AAA',native_inferred_fr3_length=1)for i in range(1,13)])
    q=pl.DataFrame({'peptide_id':['qi','ql','cross','other'],'seq':['AIGTSC','ALGTSC','AWFGAA','CAAWFG']})
    hits=pl.DataFrame([dict(chain_id=f'h{i}',peptide_id=key,start=start,end=end)for i in range(1,13)for key,start,end in [('qi',0,6),('ql',0,6),('cross',7,13),('other',5,11)]])
    blocks={BLOCKS[0]:pl.DataFrame({'seq':['AIGTSC','ALGTSC','AWFGAA','AAA','AXXX']}),BLOCKS[1]:pl.DataFrame({'seq':['AWFGAA']}),BLOCKS[2]:pl.DataFrame({'seq':['CAAWFG']})}
    native={b:[{}]*n for b,n in zip(BLOCKS,[7,5,3])}
    regions=supplied_reference_regions(panel);lengths={r['region']:r['end']-r['start']for r in regions.filter(pl.col('chain_id')=='q1').to_dicts()}
    rows=[]
    for condition in ['solution_mixture','gel_mixture']:
        for enzyme in ['trypsin','chymotrypsin','pooled']:
            rows.append(dict(condition=condition,enzyme=enzyme,reference='all',region=None,quantity='compatible_called_species',units='species',value=2.,numerator=2,denominator=None))
            for reference in ['all',*[f'ref{i}'for i in range(1,13)]]:
                for region,length in lengths.items():
                    rows.append(dict(condition=condition,enzyme=enzyme,reference=reference,region=region,quantity='residue_union_coverage',units='fraction',value=0.,numerator=0,denominator=length*(12 if reference=='all' else 1)))
    targets=pl.DataFrame(rows)
    frames,report=experimental_readouts(panel,physical,q,hits,blocks,native,targets)
    assert report['summary']['yield_counts']==[1,0,1]
    assert report['coordinates']==255 and report['shared_chymotrypsin_selection_count']==1 and not report['recovery_used_for_selection']
    cov=frames['solution_trypsin_supplied_query_coverage']
    assert cov.filter((pl.col('chain_id')=='ref1')&(pl.col('region')=='FR3'))['region_length'].item()==2
    assert frames['heavy_band_trypsin_supplied_query_coverage']['covered_residues'].sum()==0
    inv=dict(frames['solution_trypsin_query_species'].select('quantity','count').iter_rows())
    assert inv['eligible_species']==2 and inv['compatible_species']==1 and inv['unmapped_species']==1
    assert inv['short_rows']==2 and inv['noncanonical_rows']==1 and inv['native_ms2_groups']==7
    pooled=dict(frames['solution_pooled_query_species'].select('quantity','count').iter_rows())
    assert pooled['native_ms2_groups']==10 and pooled['compatible_species']==2
    residual=frames['experimental_target_residuals']
    species=residual.filter((pl.col('quantity')=='compatible_called_species')&(pl.col('enzyme')=='trypsin'))
    assert dict(species.select('condition','simulated_value').iter_rows())=={'gel_mixture':0.,'solution_mixture':1.}
    shared=residual.filter((pl.col('quantity')=='compatible_called_species')&(pl.col('enzyme')=='chymotrypsin'))
    assert shared['simulated_value'].to_list()==[1.,1.]
    # New recorded observations alter only reported residual/distance, never simulated outputs.
    changed=targets.with_columns(pl.when(pl.col('quantity')=='compatible_called_species').then(3.).otherwise(pl.col('value')).alias('value'),pl.when(pl.col('quantity')=='compatible_called_species').then(3).otherwise(pl.col('numerator')).alias('numerator'))
    _,other=experimental_readouts(panel,physical,q,hits,blocks,native,changed)
    assert other['summary']==report['summary'] and other['distance']!=report['distance']
