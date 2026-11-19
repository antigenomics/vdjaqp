"""The fresh-run preflight must reject before launching native work."""
import importlib.util
from pathlib import Path
import pytest

spec = importlib.util.spec_from_file_location('metric_driver', Path(__file__).parents[1]/'scripts/smoke_metrics.py')
module = importlib.util.module_from_spec(spec)

def test_preflight(tmp_path, monkeypatch):
    spec.loader.exec_module(module)
    monkeypatch.setattr(module.subprocess, 'run', lambda *a, **k: pytest.fail('native launched'))
    with pytest.raises(FileExistsError): module.preflight(tmp_path, 3)
    missing = tmp_path/'fresh'
    with pytest.raises(ValueError, match='positive'): module.preflight(missing, 0)
    assert not missing.exists()


def test_undefined_aggregate_stays_null():
    spec.loader.exec_module(module)
    import polars as pl
    frame = pl.DataFrame({'locus':['IGK','IGH'], 'region':['D','D'],
        'region_status':['not_applicable','available'], 'coverage':[None,0.],
        'covered_residues':[None,0], 'region_length':[None,3]})
    result = module.aggregate_coverage(frame,['locus','region'])
    light = result.filter(pl.col('locus')=='IGK').row(0,named=True)
    assert light['defined_chains']==0 and light['status']=='not_applicable'
    assert all(light[k] is None for k in ['covered_chains','residue_denominator','covered_residues','pooled_residue_coverage'])
    heavy = result.filter(pl.col('locus')=='IGH').row(0,named=True)
    assert heavy['covered_residues']==0 and heavy['residue_denominator']==3
    assert heavy['status']=='defined' and heavy['pooled_residue_coverage']==0.


def test_gene_summary_preserves_native_bins_and_ambiguity():
    spec.loader.exec_module(module)
    import polars as pl
    frame=pl.DataFrame({'locus':['IGH']*3,'region':['V','V','C'],
        'gene_call':['old*1','old*2','IGHG1'],'evidence_level':['called']*3,
        'region_status':['available']*3,'coverage':[.5,.25,1.],
        'covered_residues':[2,1,4],'region_length':[4]*3})
    metadata=pl.DataFrame({'locus':['IGH']*2,'v_call':['old*1','old*2'],
        'canonical_v_gene':['IGHV1-18',None],
        'canonical_v_match_status':['unique_gene','ambiguous_genes'],
        'canonical_v_gene_count':[1,2]})
    result=module.aggregate_gene_coverage(frame,metadata)
    assert result['represented_chains'].sum()==3
    assert result['residue_denominator'].sum()==12
    rows={r['gene_call']:r for r in result.to_dicts()}
    assert rows['old*1']['canonical_v_gene']=='IGHV1-18'
    assert rows['old*2']['canonical_v_gene'] is None
    assert rows['old*2']['canonical_v_match_status']=='ambiguous_genes'
    assert rows['IGHG1']['canonical_v_match_status']=='not_applicable'
    assert rows['IGHG1']['canonical_v_gene_count'] is None
    assert result.write_csv(separator='\t')


def test_empty_candidates_keep_physical_identity_denominator():
    import polars as pl
    from vdjaqp.bench import rank_targets
    from vdjaqp.simulation_metrics import identity_metrics
    metadata = pl.DataFrame({'chain_id':['h'],'locus':['IGH'],
        'variable_aa':['ACDEFI'],'protein_aa':['ACDEFIQQ']})
    candidates = pl.DataFrame(schema={'clone_id':pl.String})
    scores = pl.DataFrame(schema={'clone_id':pl.String,'score':pl.Float64,'confirmed_nats':pl.Float64})
    ranks = rank_targets(scores,'score',set())
    associations = pl.DataFrame(schema={'clone_id':pl.String,'chain_id':pl.String,'start':pl.Int32,'end':pl.Int32})
    population, called, summary = identity_metrics(metadata,candidates,ranks,associations)
    assert population.height==1 and called.height==0
    row = summary.filter(pl.col('quantity')=='recovered_classes_at20').row(0,named=True)
    assert row['count']==0 and row['denominator']==1


@pytest.mark.parametrize('flags', [
    ['--families','0'], ['--families','3'], ['--families','10000000000'],
    ['--pairs-per-chain','100'], ['--mutation-rate','.01'],
    ['--protein-concentration','equal'], ['--seed','123'],
    ['--families','40','--pairs-per-chain','0'],
    ['--families','40','--mutation-rate','nan'],
    ['--families','40','--mutation-rate','inf'],
    ['--families','40','--mutation-rate','1.1'],
    ['--families','40','--seed','-1'],
    ['--gradient-seconds','nan'], ['--gradient-seconds','0'],
    ['--ms2-per-bin','0'],
])
def test_scenario_preflight_before_output_or_native(tmp_path, monkeypatch, flags):
    spec.loader.exec_module(module)
    output = tmp_path/'absent'
    monkeypatch.setattr(module.sys, 'argv', ['smoke_metrics.py','--output',str(output),*flags])
    monkeypatch.setattr(module.subprocess, 'run', lambda *a, **k: pytest.fail('native launched'))
    with pytest.raises(SystemExit) as error:
        module.main()
    assert error.value.code == 2
    assert not output.exists()


def test_population_defaults_and_legacy_route(tmp_path):
    spec.loader.exec_module(module)
    flags = ['--output',str(tmp_path/'fresh')]
    legacy = module.parse_args(flags)
    assert legacy.families is None and legacy.pairs_per_chain is None
    assert legacy.gradient_seconds == 120 and legacy.ms2_per_bin == 5
    args = module.parse_args([*flags,'--families','40','--purpose','resource-profile',
        '--gradient-seconds','1200','--ms2-per-bin','10'])
    assert (args.families,args.pairs_per_chain,args.mutation_rate,args.protein_concentration,args.seed) == (40,100,.01,'equal',20261110)
    assert args.purpose == 'resource-profile' and args.gradient_seconds == 1200 and args.ms2_per_bin == 10


@pytest.mark.parametrize('explicit_raw', [False, True])
@pytest.mark.parametrize('heavy_gene', ['IGHG1','IGHM'])
def test_population_generation_keeps_output_absent_and_only_bootstraps_raw(tmp_path, monkeypatch, explicit_raw, heavy_gene):
    import sys
    from types import SimpleNamespace
    spec.loader.exec_module(module)
    output, raw = tmp_path/'fresh', tmp_path/'raw'
    flags = ['--output',str(output),'--families','40']
    additional = tmp_path/'additional'
    if heavy_gene != 'IGHG1':
        additional.mkdir()
        flags += ['--heavy-constant-gene',heavy_gene,'--additional-reference-dir',str(additional)]
    if explicit_raw:
        raw.mkdir()
        flags += ['--reference-dir',str(raw)]
    prepared, launched = [], []
    def prepare(repo, root, **kwargs):
        assert not explicit_raw and root != output
        prepared.append(root)
        (root/'references').mkdir()
    monkeypatch.setitem(sys.modules,'smoke_all',SimpleNamespace(prepare_references=prepare))
    def run(label, interpreter, script, *arguments):
        assert not output.exists()
        assert label == 'generation' and script == 'generate_population.py'
        pairs = dict(zip(arguments[::2], arguments[1::2]))
        reference = Path(pairs['--reference-dir'])
        assert reference.exists()
        assert pairs['--output'] == output and pairs['--families'] == 40
        assert pairs['--heavy-constant-gene'] == heavy_gene
        if heavy_gene != 'IGHG1':
            assert pairs['--additional-reference-dir'] == additional.resolve()
        else:
            assert '--additional-reference-dir' not in pairs
        launched.append(reference)
    module.generate_inputs(module.parse_args(flags),tmp_path,output,run)
    assert len(launched) == 1 and not output.exists()
    if explicit_raw:
        assert launched[0] == raw and raw.exists() and not prepared
    else:
        assert len(prepared) == 1 and not prepared[0].exists()


def test_heavy_constant_options_are_explicit_and_population_only(tmp_path):
    spec.loader.exec_module(module)
    raw=tmp_path/'additional';raw.mkdir()
    base=['--output',str(tmp_path/'fresh')]
    for extra in [ ['--heavy-constant-gene','IGHM'],
                   ['--families','2','--heavy-constant-gene','IGHM'],
                   ['--families','2','--additional-reference-dir',str(raw)] ]:
        with pytest.raises(SystemExit):module.parse_args(base+extra)
    args=module.parse_args(base+['--families','2','--heavy-constant-gene','IGHM',
                                '--additional-reference-dir',str(raw)])
    assert args.heavy_constant_gene=='IGHM' and args.additional_reference_dir==raw
    assert module.parse_args(base+['--families','2']).heavy_constant_gene=='IGHG1'


def test_search_constants_are_declared_reference_inputs(tmp_path):
    spec.loader.exec_module(module)
    with pytest.raises(FileNotFoundError):module.constant_search_fasta(tmp_path,population=True)
    text='>constant_selected\nACDEFG\n'
    (tmp_path/'constant_reference_proteins.fasta').write_text(text)
    (tmp_path/'truth.json').write_text('{"protein_aa":"NEVER_USE_GENERATED_VARIABLE_TRUTH"}')
    assert module.constant_search_fasta(tmp_path,population=True)==text
    raw=tmp_path/'references';raw.mkdir()
    accessions=('P01857-1','P01834','P0DOY2')
    for accession in accessions:(raw/(accession+'.fasta')).write_text('>source\nACDEFG\n')
    expected=''.join('>constant_'+a+'\nACDEFG\n' for a in accessions)
    assert module.constant_search_fasta(tmp_path,population=False)==expected


def test_mixture_route_and_component_denominators(tmp_path, monkeypatch):
    import sys
    from types import SimpleNamespace
    import polars as pl
    spec.loader.exec_module(module)
    raw=tmp_path/'raw';raw.mkdir()
    additional=tmp_path/'additional';additional.mkdir()
    flags=['--output',str(tmp_path/'fresh'),'--total-antibodies','1200',
           '--memory-fraction','.5','--reference-dir',str(raw),
           '--additional-reference-dir',str(additional)]
    args=module.parse_args(flags)
    assert args.families is None and args.total_antibodies==1200 and args.memory_fraction==.5
    for extra in [['--families','200'],['--heavy-constant-gene','IGHM']]:
        with pytest.raises(SystemExit):module.parse_args(flags+extra)
    for invalid in [['--total-antibodies','1200'],['--total-antibodies','14','--memory-fraction','.5']]:
        with pytest.raises(SystemExit):module.parse_args(['--output',str(tmp_path/'fresh'),*invalid])
    monkeypatch.setitem(sys.modules,'smoke_all',SimpleNamespace(prepare_references=lambda *a,**kw:pytest.fail('unexpected download')))
    commands=[]
    def run(label, interpreter, script, *arguments):
        command=dict(zip(arguments[::2],arguments[1::2]))
        assert '--families' not in command
        assert command['--total-antibodies']==1200 and command['--memory-fraction']==.5
        assert command['--additional-reference-dir']==additional.resolve()
        commands.append(command)
    module.generate_inputs(args,tmp_path,args.output,run)
    assert len(commands)==1 and not args.output.exists()
    metadata=pl.DataFrame({'chain_id':['h1','h2','h3'], 'locus':['IGH']*3,
        'component':['naive','memory','memory'], 'antibody_id':['a1','a2','a3'],
        'family_id':['f1','f2','f2']})
    counts=module.component_counts(metadata)
    assert counts['observed_chain_count'].sum()==3
    memory=counts.filter(pl.col('component')=='memory').row(0,named=True)
    assert memory['observed_antibody_count']==2 and memory['observed_family_count']==1
    identity=metadata.with_columns(pl.Series('variable_mass_class_id',['same','same','other']),
        pl.Series('full_variable_compatible_candidate_count',[1,1,0]),
        pl.Series('class_recovered_at20',[False,False,False]))
    recovery=module.component_recovery(identity)
    assert recovery['observed_heavy_rows'].sum()==3
    assert recovery['recovered_classes_at20'].sum()==0
    memory=recovery.filter(pl.col('component')=='memory').row(0,named=True)
    assert memory['variable_mass_classes']==2 and memory['observable_classes']==1
    assert memory['class_recovery_at20']==0 and memory['conditional_class_recovery_at20']==0
    assert counts.write_csv(separator='\t') and recovery.write_csv(separator='\t')


@pytest.mark.parametrize('geometry', [None, 'both', 'amplicon', 'rnaseq'])
@pytest.mark.parametrize('timeout', [None, 300., 900.])
def test_geometry_forwarded_only_when_explicit(tmp_path, monkeypatch, geometry, timeout):
    import json
    spec.loader.exec_module(module)
    root = tmp_path/'fresh'
    flags = ['smoke_metrics.py', '--output', str(root), '--cpus', '3']
    if geometry is not None:
        flags += ['--geometry', geometry]
    if timeout is not None: flags += ['--stage-timeout', str(timeout)]
    monkeypatch.setattr(module.sys, 'argv', flags)
    from test_simulation_truth import fixture
    import arda.annotate.contig
    import vdjaqp.simulation_truth
    chain, annotation = fixture()
    chain.update(rna_pairs=1,secretion_factor=1.,family_id='family0',antibody_id='antibody0')
    events = []
    modes=module.selected_modes(geometry or 'both')
    def annotate(queries, **kwargs):
        assert events == ['reads']+[mode+'_repertoire' for mode in modes]
        events.append('annotation')
        assert queries == [(chain['chain_id'], chain['variable_nt'])]
        # --cpus 3 is clamped to the host's available cores (2 on CI runners).
        assert kwargs == dict(seqtype='nt', threads=vdjaqp.resolve_cpu_budget(3))
        return [annotation]
    original = vdjaqp.simulation_truth.regions_from_annotations
    def validate(*args):
        events.append('validated')
        return original(*args)
    monkeypatch.setattr(arda.annotate.contig, 'reannotate_contigs', annotate)
    monkeypatch.setattr(vdjaqp.simulation_truth, 'regions_from_annotations', validate)
    def inputs(*args):
        root.mkdir()
        (root/'read_counts.tsv').write_text('chain\t2\n')
        (root/'truth.json').write_text(json.dumps({'chains':[chain]}))
    monkeypatch.setattr(module, 'generate_inputs', inputs)
    commands = []
    class FinishedObservations(Exception): pass
    def native(command, **kwargs):
        commands.append(command)
        assert kwargs['timeout'] == (timeout or 300.)
        assert ('--stage-timeout' in command) == (timeout is not None)
        if timeout is not None: assert command[command.index('--stage-timeout')+1] == str(timeout)
        if command[1].endswith('smoke_reads.py'):
            assert '--geometry' in command if geometry is not None else '--geometry' not in command
            if geometry is not None: assert command[command.index('--geometry')+1] == geometry
            events.append('reads')
            (root/'read_run_manifest.json').write_text(json.dumps([
                {'scenario':mode,'pairs':1} for mode in modes]))
        else:
            assert command[1].endswith('smoke_repertoire.py')
            assert command[2] in modes
            assert '--geometry' not in command
            events.append(command[2]+'_repertoire')
    def load(path, mode):
        assert path==root/(mode+'_observed')
        assert events == ['reads']+[m+'_repertoire' for m in modes]+['annotation','validated']
        raise FinishedObservations
    monkeypatch.setattr(module, 'load_observed', load)
    monkeypatch.setattr(module.subprocess, 'run', native)
    with pytest.raises(FinishedObservations):
        module.main()
    assert len(commands) == 1+len(modes)
    assert module.parse_args(['--output',str(tmp_path/'other')]).geometry is None


def test_explicit_simulator_and_engineering_timeout_preflight(tmp_path):
    spec.loader.exec_module(module)
    flags=['--output',str(tmp_path/'fresh')]
    assert module.parse_args(flags).stage_timeout==300.
    binary=tmp_path/'MSSimulator-fixed';binary.write_text('#!/bin/sh\nexit 0\n')
    with pytest.raises(SystemExit):module.parse_args(flags+['--mssimulator',str(binary)])
    binary.chmod(0o755)
    args=module.parse_args(flags+['--mssimulator',str(binary),'--stage-timeout','900'])
    assert args.mssimulator==binary.resolve() and args.stage_timeout==900
    for timeout in ('0','-1','nan','inf'):
        with pytest.raises(SystemExit):module.parse_args(flags+['--stage-timeout',timeout])


def test_generation_timeout_and_reads_before_truth(tmp_path, monkeypatch):
    import json
    import arda.annotate.contig
    from test_simulation_truth import fixture
    spec.loader.exec_module(module)
    root = tmp_path/'fresh'
    monkeypatch.setattr(module.sys,'argv',['smoke_metrics.py','--output',str(root),
        '--cpus','3','--stage-timeout','900'])
    chain, annotation = fixture()
    chain['rna_pairs'] = 1
    def inputs(args, repo, output, run):
        run('generation', module.sys.executable, 'smoke_generate.py', '--output', root)
    monkeypatch.setattr(module,'generate_inputs',inputs)
    commands=[]
    class ReachedReads(Exception): pass
    def native(command, **kwargs):
        assert kwargs['timeout']==900.
        commands.append(command)
        if command[1].endswith('smoke_reads.py'):
            assert command[command.index('--stage-timeout')+1]=='900.0'
            raise ReachedReads
        assert command[1].endswith('smoke_generate.py')
        assert '--stage-timeout' not in command
        root.mkdir()
        (root/'read_counts.tsv').write_text('chain\t2\n')
        (root/'truth.json').write_text(json.dumps({'chains':[chain]}))
    monkeypatch.setattr(module.subprocess,'run',native)
    monkeypatch.setattr(arda.annotate.contig,'reannotate_contigs',lambda *a,**k:pytest.fail('truth annotation ran before RNA'))
    with pytest.raises(ReachedReads):
        module.main()
    assert len(commands)==2
    assert not (root/'truth_annotations.json').exists()


@pytest.mark.parametrize('geometry',['amplicon','rnaseq'])
def test_native_export_receipt_checks(tmp_path,geometry):
    import hashlib,json
    import polars as pl
    spec.loader.exec_module(module)
    source=tmp_path/'native.airr.tsv';source.write_text('constructed native source\n')
    receipt=dict(status='complete',source_format='arda_fastq_native',geometry=geometry,counts={},artifact_sha256={},
        source_input_sha256={str(source):hashlib.sha256(source.read_bytes()).hexdigest(),
                             str(tmp_path/'unavailable_remote.tsv'):'0'*64})
    for name in ('catalog','candidates','search_entries'):
        file=tmp_path/(name+'.parquet')
        pl.DataFrame({'background_member':[True],'native_query':['unchanged']}).write_parquet(file)
        receipt['counts'][name]=1
        receipt['artifact_sha256'][file.name]=hashlib.sha256(file.read_bytes()).hexdigest()
    (tmp_path/'receipt.json').write_text(json.dumps(receipt))
    observed=module.load_observed(tmp_path,geometry)
    assert observed['candidates']['native_query'].to_list()==['unchanged']
    assert observed['source_verification']==dict(verified_local_sources=[str(source)],
        unavailable_sources_not_verified=[str(tmp_path/'unavailable_remote.tsv')])
    with pytest.raises(ValueError,match='geometry'):
        module.load_observed(tmp_path,'wrong')
    for key in ('hash','count','background','status','source'):
        changed=json.loads(json.dumps(receipt))
        if key=='hash': changed['artifact_sha256']['candidates.parquet']='invalid'
        if key=='count': changed['counts']['candidates']=2
        if key=='status': changed['status']='incomplete'
        if key=='source': changed['source_input_sha256'][str(source)]='0'*64
        if key=='background':
            file=tmp_path/'candidates.parquet'
            pl.DataFrame({'background_member':[False]}).write_parquet(file)
            changed['artifact_sha256'][file.name]=hashlib.sha256(file.read_bytes()).hexdigest()
        (tmp_path/'receipt.json').write_text(json.dumps(changed))
        with pytest.raises(ValueError,match={'hash':'artifact hash','count':'count mismatch',
            'background':'membership','status':'must be complete','source':'source hash'}[key]):
            module.load_observed(tmp_path,geometry)
