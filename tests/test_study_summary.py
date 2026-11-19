"""Final-table summaries: exact inventory/provenance and no scientific execution."""
import gzip
import importlib.util
import json
from pathlib import Path
import polars as pl
import pytest

spec=importlib.util.spec_from_file_location('study_summary',Path(__file__).parents[1]/'scripts/summarize_simulation_study.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def mock_bundle(root, case='M0', mode='amplicon'):
    root.mkdir(parents=True)
    fraction,mutation,concentration,capacity=m.CASES[case]
    naive=1200-int(1200*fraction)
    rows=[]
    for locus,n in [('IGH',1200),('IGK',600),('IGL',600)]:
        for i in range(n):
            index=i if locus=='IGH' else i*2
            rows.append(dict(chain_id=f'{locus}{i}',antibody_id=f'a{index}',family_id=f'f{index}',locus=locus,
                component='naive' if index<naive else 'memory',rna_pairs=100,serum_weight=2.,
                variable_aa='ACDEFI',protein_aa='ACDEFIPP',generator_junction_nt='TGT',generator_junction_aa='CAAW',
                annotated_junction_nt='TGT',annotated_junction_aa='CAAW',v_call='nativeV',j_call='nativeJ',
                constant_nt_gene='IGHM' if index<naive else 'IGHG1',canonical_v_gene=None,
                canonical_v_match_status='no_exact_match'))
    pl.DataFrame(rows,schema_overrides={'canonical_v_gene':pl.String}).write_parquet(root/'repertoire_metrics.parquet')
    candidate=pl.DataFrame(dict(clone_id=['c'],seq=['ACDEFI'],nt_seq=['GCT'],cdr3_seq=['CAAW'],
        region_ids=[[0]*6],v_name=['observedV'],j_name=['observedJ'],c_name=[None]),schema_overrides={'c_name':pl.String})
    candidate.write_parquet(root/(mode+'_candidate_identity.parquet'))
    pl.DataFrame(dict(sequence_id=['c'],locus=['IGH'],status=['usable_heavy_candidate'])).write_parquet(root/'contig_catalog.parquet')
    annotations=[dict(sequence_id='c',locus='IGH',v_call='observedV',j_call='observedJ',c_call=None,junction_aa='CAAW')]
    assembly=root/'assembly';assembly.mkdir()
    raw=json.dumps(annotations).encode();(assembly/'arda_contigs.json.gz').write_bytes(gzip.compress(raw,mtime=0))
    dummy=pl.DataFrame(dict(quantity=['existing_scalar'],count=[1]))
    for name in [*(''+mode+'_'+x for x in m.EXISTING),'truth_regions','component_counts']:
        dummy.write_parquet(root/(name+'.parquet'))
    for enzyme,seqs in [('Trypsin',['ACDEFI','ACDEFL']),('Chymotrypsin',['ACDEFL','CAAAAA'])]:
        pl.DataFrame(dict(seq=seqs,enzyme=[enzyme]*2,correct=[True,False])).write_parquet(root/(mode+'_'+enzyme+'_psm_correctness.parquet'))
        dummy.write_parquet(root/(mode+'_'+enzyme+'_readout_metrics.parquet'))
    scientific={name:'a'*64 for name in m.SCIENTIFIC_ENTRYPOINTS}
    monitor='57d57da451927c7f695fdb2f1d0387ac6a7fe014427bdb9281b36c07919b040a'
    (root/'entrypoints.sha256').write_text(''.join(sha+'  '+name+'\n' for name,sha in sorted(scientific.items()))
        +monitor+'  scripts/profile_build.py\n')
    resource=dict(failure=None,exit_code=0,cpu_budget=16,live_thread_limit=64,peak_live_threads=10,
        peak_sampled_summed_rss_bytes=100,wall_seconds=1,peak_sampled_project_allocated_bytes=100)
    (root/'resource.json').write_text(json.dumps(resource))
    params=dict(geometry=mode,memory_fraction=fraction,mutation_rate=mutation,protein_concentration=concentration,
        ms2_per_bin=capacity,seed=123,pairs_per_chain=100,gradient_seconds=1200,purpose='mechanistic',
        reconstruction_claim=False,antibodies=1200,chains=2400,total_antibodies=1200,actual_read_pairs={mode:240000})
    files=[p for p in root.rglob('*') if p.is_file()]
    hashes={str(p.relative_to(root)):m.digest(p) for p in files if p.name!='arda_contigs.json.gz'}
    hashes['assembly/arda_contigs.json']=__import__('hashlib').sha256(raw).hexdigest()
    report=dict(scenario_inputs=params,scenarios={mode:{}},artifact_sha256=hashes,
        discarded_spectra={f'{mode}_{e}/spectra.mzML':{} for e in ('Trypsin','Chymotrypsin')},
        software_sha256=scientific)
    (root/'metric_integration.json').write_text(json.dumps(report))
    records=[]
    for p in [*files,root/'metric_integration.json']:
        name=str(p.relative_to(root));record=dict(file=name,bytes=p.stat().st_size,sha256=m.digest(p),source_path=name)
        if name.endswith('.gz'):
            record.update(source_path='assembly/arda_contigs.json',source_sha256=hashes['assembly/arda_contigs.json'],
                source_bytes=len(raw),compression='gzip; mtime=0; empty header filename')
        records.append(record)
    (root/'retained_table_manifest.json').write_text(json.dumps(dict(artifacts=records)))
    return records,report


def test_exact_inventory_hashes_and_mock_final_summaries(tmp_path):
    inputs=tmp_path/'inputs';inputs.mkdir()
    with pytest.raises(ValueError,match='twelve'):m.inventory(inputs)
    for case in m.CASES:
        for mode in m.GEOMETRIES: mock_bundle(inputs/f'{case}-{mode}',case,mode)
    assert len(m.inventory(inputs))==12
    root=inputs/'M0-amplicon'
    records,report,monitor=m.validate_bundle(root,'M0','amplicon')
    assert monitor['scripts/profile_build.py'] in m.MONITOR_PINS
    assert 'scripts/profile_build.py' not in report['software_sha256']
    result=m.summarize_tables(root,'M0','amplicon',records)
    assert result['generated_antibodies'].height==1200
    assert result['generated_antibodies']['protein_antibody_fraction'].sum()==pytest.approx(1.)
    assert result['generated_antibodies']['rna_antibody_fraction'].sum()==pytest.approx(1.)
    assert result['generated_chains']['generator_cdr3_length_aa'].unique().to_list()==[2]
    counts=dict(zip(result['peptide_species']['enzyme'],result['peptide_species']['count']))
    assert counts=={'Trypsin':1,'Chymotrypsin':2,'pooled':2}
    assert result['pooled_PSM_accuracy']['numerator'].item()==2
    assert result['pooled_PSM_accuracy']['denominator'].item()==4
    genes=result['gene_usage']
    assert set(genes['gene_namespace'])=={'native_GenAIRR_OGRDB','native_ARDA_calls'}
    assert genes.filter((pl.col('population')=='observed_raw_contig_records')&(pl.col('segment')=='C'))['missing_call'].item()
    for frame in result.values():
        assert not any(t.is_nested() for t in frame.schema.values())
        frame.write_csv(separator='\t')
    # Exercise the complete 12-folder CLI using only synthetic final-table fixtures.
    output=tmp_path/'summary'
    m.main(['--input-root',str(inputs),'--output',str(output),'--cpus','1'])
    provenance=json.loads((output/'provenance.json').read_text())
    assert len(provenance['cases'])==12
    assert provenance['cases'][0]['monitoring_source_sha256']==monitor
    assert result['identity_metrics']['population'].unique().to_list()==['producer_defined_heavy_identity_readouts']
    assert result['identity_population']['population'].unique().to_list()==['physical_heavy_records']
    repeated=tmp_path/'summary-repeat'
    m.main(['--input-root',str(inputs),'--output',str(repeated),'--cpus','1'])
    assert {p.name:p.read_bytes() for p in output.glob('*.tsv.gz')} == {p.name:p.read_bytes() for p in repeated.glob('*.tsv.gz')}
    assert (output/'provenance.json').read_bytes()==(repeated/'provenance.json').read_bytes()
    with pytest.raises(SystemExit):m.main(['--input-root',str(inputs),'--output',str(output)])
    (inputs/'extra').mkdir()
    with pytest.raises(ValueError,match='twelve'):m.inventory(inputs)
    (inputs/'extra').rmdir()
    path=root/'metric_integration.json';path.write_text(path.read_text()+' ')
    with pytest.raises(ValueError,match='hash/size'):m.validate_bundle(root,'M0','amplicon')


def test_wrong_parameters_and_manifest_escape_rejected(tmp_path):
    root=tmp_path/'case';mock_bundle(root)
    path=root/'metric_integration.json';r=json.loads(path.read_text());r['scenario_inputs']['seed']=999;path.write_text(json.dumps(r))
    manifest=root/'retained_table_manifest.json';rows=json.loads(manifest.read_text())
    rec=next(a for a in rows['artifacts'] if a['file']=='metric_integration.json')
    rec.update(bytes=path.stat().st_size,sha256=m.digest(path));manifest.write_text(json.dumps(rows))
    with pytest.raises(ValueError,match='parameters'):m.validate_bundle(root,'M0','amplicon')
    rec['file']='../escape.json';manifest.write_text(json.dumps(rows))
    with pytest.raises(ValueError):m.validate_bundle(root,'M0','amplicon')


def test_monitor_exception_is_explicit_and_scientific_sources_stay_required(tmp_path):
    for variant in ('unknown_monitor','missing_scientific','unreported_source','new_monitor'):
        root=tmp_path/variant;mock_bundle(root)
        path=root/'entrypoints.sha256';text=path.read_text()
        if variant=='unknown_monitor':
            text=text.replace('57d57da451927c7f695fdb2f1d0387ac6a7fe014427bdb9281b36c07919b040a','f'*64)
        elif variant=='missing_scientific':
            text='\n'.join(line for line in text.splitlines() if not line.endswith('scripts/smoke_reads.py'))+'\n'
        elif variant=='unreported_source':
            text+='a'*64+'  scripts/unreported.py\n'
        else:
            text=text.replace('57d57da451927c7f695fdb2f1d0387ac6a7fe014427bdb9281b36c07919b040a',
                '39895acf682ef491f0d01ff178a83fa73ed3dbdef458ec229f1e63300d147acf')
        path.write_text(text)
        report_path=root/'metric_integration.json';report=json.loads(report_path.read_text())
        report['artifact_sha256']['entrypoints.sha256']=m.digest(path)
        report_path.write_text(json.dumps(report))
        manifest_path=root/'retained_table_manifest.json';manifest=json.loads(manifest_path.read_text())
        for r in manifest['artifacts']:
            if r['file'] in ('entrypoints.sha256','metric_integration.json'):
                file=root/r['file'];r.update(bytes=file.stat().st_size,sha256=m.digest(file))
        manifest_path.write_text(json.dumps(manifest))
        if variant=='new_monitor':
            _,_,monitor=m.validate_bundle(root,'M0','amplicon')
            assert monitor['scripts/profile_build.py'].startswith('39895acf')
        else:
            with pytest.raises(ValueError):m.validate_bundle(root,'M0','amplicon')


def transcript_bundle(root, case='M0'):
    """Only constructed source/final-table receipts, with the new packet layout."""
    mock_bundle(root,case,'rnaseq')
    report=json.loads((root/'metric_integration.json').read_text())
    params=report['scenario_inputs'];params.update(rna_model='leader_imputed_v1',stage_timeout_seconds=900.)
    resources=root/'resources';resources.mkdir()
    resource=json.loads((root/'resource.json').read_text())
    resource['command']=['/mock/source-9aa3efe/.venv/bin/python','/mock/source-9aa3efe/scripts/smoke_metrics.py']
    (resources/'resource.json').write_text(json.dumps(resource))
    (resources/'validation.json').write_text(json.dumps(dict(case=case,status='verified_final_deliverables',
        scenario_inputs=params,resource=resource,verified_source_count=len(report['software_sha256']))))
    names=set(report['software_sha256'])|{'scripts/profile_build.py'}
    (resources/'verified_source_files.txt').write_text(''.join('source-9aa3efe/'+n+': OK\n' for n in sorted(names)))
    (resources/'cli_help.txt').write_text('fixture')
    (resources/'job.log').write_text('fixture')
    (root/'resource.json').unlink();(root/'entrypoints.sha256').unlink()
    report['artifact_sha256'].pop('resource.json');report['artifact_sha256'].pop('entrypoints.sha256')
    report['stages']=[dict(stage=stage) for stage in ['generation','reads','assembly','evaluation_truth_annotation',
        'observed_transcript_boundaries','rnaseq_Trypsin','rnaseq_Chymotrypsin']]
    (root/'read_run_manifest.json').write_text(json.dumps([dict(scenario='rnaseq',pairs=240000,model='MiSeq',cpus=16,
        command=['iss','generate','--seed','2468','--fragment-length','400','--fragment-length-sd','40'])]))
    truth=dict(rna_model='leader_imputed_v1',transcript_model={'model_source_sha256':m.TRANSCRIPT_MODEL_SHA256})
    (root/'truth.json.gz').write_bytes(gzip.compress(json.dumps(truth).encode(),mtime=0))
    pl.DataFrame(dict(sequence_id=['c'],status=['resolved'])).write_parquet(root/'observed_transcript_boundaries.parquet')
    for name in ['read_run_manifest.json','observed_transcript_boundaries.parquet']:
        report['artifact_sha256'][name]=m.digest(root/name)
    report['artifact_sha256']['truth.json']=__import__('hashlib').sha256(json.dumps(truth).encode()).hexdigest()
    (root/'metric_integration.json').write_text(json.dumps(report))
    refresh_manifest(root,case)
    (root.parent/(case+'_mature_payload_invariant.json')).write_text(json.dumps(dict(case=case,status='pass',
        n_physical_chains_old=2400,n_physical_chains_new=2400,old_physical_payload_sha256='a'*64,
        new_physical_payload_sha256='a'*64,new_truth_sha256=report['artifact_sha256']['truth.json'],
        RNA_prefix_only=True,latent_ancestors_equal=True,new_rna_model='leader_imputed_v1',first_mismatch=None,
        top_generation_fields_equal={'seed':True})))
    return report


def refresh_manifest(root,case):
    records=[]
    for p in sorted(root.rglob('*')):
        if not p.is_file() or p.name=='retained_table_manifest.json':continue
        name=str(p.relative_to(root));rec=dict(file=name,bytes=p.stat().st_size,sha256=m.digest(p),source_path=name)
        if name.endswith('.gz'):
            raw=gzip.decompress(p.read_bytes())
            rec.update(source_path=name[:-3],source_bytes=len(raw),source_sha256=__import__('hashlib').sha256(raw).hexdigest(),compression='gzip')
        records.append(rec)
    (root/'retained_table_manifest.json').write_text(json.dumps(dict(case=case,model='leader_imputed_v1',artifacts=records)))


def test_explicit_transcript_mode_six_inventory_and_same_readout_engine(tmp_path,monkeypatch):
    root=tmp_path/'input';root.mkdir()
    with pytest.raises(ValueError,match='six'):m.inventory(root,'leader_imputed_v1')
    for case in m.CASES:report=transcript_bundle(root/(case+'-rnaseq'),case)
    monkeypatch.setattr(m,'TRANSCRIPT_SOFTWARE_SHA256',m.software_digest(report['software_sha256']))
    assert len(m.inventory(root,'leader_imputed_v1'))==6
    with pytest.raises(ValueError,match='twelve'):m.inventory(root)
    output=tmp_path/'out'
    m.main(['--input-root',str(root),'--output',str(output),'--study-model','leader_imputed_v1','--cpus','1'])
    assert len(list(output.glob('*.tsv.gz')))==16
    provenance=json.loads((output/'provenance.json').read_text())
    assert provenance['study_model']=='leader_imputed_v1' and len(provenance['cases'])==6
    assert {x['producer_commit'] for x in provenance['cases']}=={'9aa3efe'}
    for p in output.glob('*.tsv.gz'):
        frame=pl.read_csv(p,separator='\t',null_values='null')
        assert frame['rna_model'].unique().to_list()==['leader_imputed_v1']
    (root/'M0-amplicon').mkdir()
    with pytest.raises(ValueError,match='six'):m.inventory(root,'leader_imputed_v1')


@pytest.mark.parametrize('defect',['model','source','read_seed','resource','truth_model','stage','payload'])
def test_transcript_payload_mismatch_fails_before_export(tmp_path,monkeypatch,defect):
    root=tmp_path/'M0-rnaseq';report=transcript_bundle(root)
    monkeypatch.setattr(m,'TRANSCRIPT_SOFTWARE_SHA256',m.software_digest(report['software_sha256']))
    if defect=='model':report['scenario_inputs']['rna_model']='coding_only'
    elif defect=='source':report['software_sha256']['scripts/smoke_reads.py']='f'*64
    elif defect=='read_seed':
        p=root/'read_run_manifest.json';p.write_text(p.read_text().replace('2468','2469'))
        report['artifact_sha256'][p.name]=m.digest(p)
    elif defect=='resource':
        p=root/'resources/resource.json';r=json.loads(p.read_text());r['exit_code']=1;p.write_text(json.dumps(r))
    elif defect=='truth_model':
        p=root/'truth.json.gz';raw=gzip.decompress(p.read_bytes()).replace(b'leader_imputed_v1',b'coding_only')
        p.write_bytes(gzip.compress(raw,mtime=0));report['artifact_sha256']['truth.json']=__import__('hashlib').sha256(raw).hexdigest()
    elif defect=='stage':report['stages']=report['stages'][1:]
    else:
        p=root.parent/'M0_mature_payload_invariant.json';v=json.loads(p.read_text());v['new_physical_payload_sha256']='b'*64;p.write_text(json.dumps(v))
    (root/'metric_integration.json').write_text(json.dumps(report));refresh_manifest(root,'M0')
    with pytest.raises(ValueError):m.validate_bundle(root,'M0','rnaseq','leader_imputed_v1')


def test_original_default_exports_match_frozen_prechange_engine(tmp_path,monkeypatch):
    # Exact source from 9aa3efe, not computed rows; works in source-only deployments.
    # Only its three package import paths changed in the igseq -> vdjaqp rename (18 Nov 2026).
    import hashlib
    prior=Path(__file__).parent/'fixtures/prechange_study_summary.py'
    assert hashlib.sha256(prior.read_bytes()).hexdigest()=='ba0e918a7c3520864f8c8935432a48c667ffc4d12ac828f1dcd59ced59d899c4'
    spec=importlib.util.spec_from_file_location('previous_study_summary',prior)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    root=tmp_path/'inputs';root.mkdir()
    for case in m.CASES:
        for geometry in m.GEOMETRIES:mock_bundle(root/(case+'-'+geometry),case,geometry)
    before=tmp_path/'before';after=tmp_path/'after'
    old.main(['--input-root',str(root),'--output',str(before),'--cpus','1'])
    m.main(['--input-root',str(root),'--output',str(after),'--cpus','1'])
    assert {p.name:p.read_bytes() for p in before.glob('*.tsv.gz')}=={p.name:p.read_bytes() for p in after.glob('*.tsv.gz')}
    a=json.loads((before/'provenance.json').read_text());b=json.loads((after/'provenance.json').read_text())
    # The producing summary script changed; all input/row/counting semantics must not.
    assert a.pop('summary_script_sha256')!=b.pop('summary_script_sha256')
    assert a==b


def test_coding_default_refuses_transcript_model_even_with_old_packet_layout(tmp_path):
    root=tmp_path/'old-layout';mock_bundle(root)
    p=root/'metric_integration.json';report=json.loads(p.read_text());report['scenario_inputs']['rna_model']='leader_imputed_v1'
    p.write_text(json.dumps(report))
    manifest=root/'retained_table_manifest.json';rows=json.loads(manifest.read_text())
    record=next(r for r in rows['artifacts'] if r['file']=='metric_integration.json')
    record.update(bytes=p.stat().st_size,sha256=m.digest(p));manifest.write_text(json.dumps(rows))
    with pytest.raises(ValueError,match='explicit --study-model'):m.validate_bundle(root,'M0','amplicon')
