"""Inspect requested tiny tool-run deliverables; never regenerate fixtures implicitly."""
import gzip
import os
import json
from pathlib import Path
import pytest

REPO = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get('VDJAQP_SMOKE_ROOT', REPO / 'tmp' / 'tool-smoke-complete'))
fixture_required = pytest.mark.skipif(not (ROOT/'truth.json').exists(), reason='run documented tool integration fixture explicitly')


@fixture_required
def test_paired_family_truth_and_constants():
    chains=json.loads((ROOT/'truth.json').read_text())['chains']
    assert len(chains)==6
    for locus in ('IGH','IGK'):
        rows=[r for r in chains if r['locus']==locus]
        assert len(rows)==3 and rows[0]['parent_id'] is None
        assert all(r['parent_id']==rows[0]['chain_id'] for r in rows[1:])
        assert all(r['variable_nt']!=rows[0]['variable_nt'] for r in rows[1:])
        expected_constant_length = 330 if locus == 'IGH' else 107
        assert all(len(r['protein_aa'])-len(r['variable_aa']) == expected_constant_length and '*' not in r['protein_aa'] for r in rows)


@fixture_required
def test_two_published_read_geometries_and_realized_pair_counts():
    runs=json.loads((ROOT/'read_run_manifest.json').read_text())
    assert [r['scenario'] for r in runs]==['amplicon','rnaseq']
    assert 'amplicon' in runs[0]['command'] and 'metagenomics' in runs[1]['command']
    for run in runs:
        assert run['pairs']==4000
        for suffix in ('R1','R2'):
            with gzip.open(ROOT/f'{run["scenario"]}_{suffix}.fastq.gz','rt') as f:
                assert sum(1 for _ in f)==16000


@fixture_required
def test_nt_constant_splice_agrees_with_full_protein():
    import subprocess
    interpreter=REPO/'.tool-envs/generation/bin/python'
    if not interpreter.exists():
        pytest.skip('published-tool translation runtime not installed')
    code = """
import json, sys
from Bio.Seq import translate
from pathlib import Path
for row in json.loads(Path(sys.argv[1]).read_text())['chains']:
    nt = row['rnaseq_nt']
    aa = str(translate(nt[:len(nt)//3*3]))
    assert aa == row['protein_aa']+'*', row['chain_id']
"""
    subprocess.run([str(interpreter),'-c',code,str(ROOT/'truth.json')],check=True)


def test_generate_refuses_existing_outputs_before_loading_references(tmp_path):
    import subprocess
    interpreter=REPO/'.tool-envs/generation/bin/python'
    if not interpreter.exists():
        pytest.skip('generation runtime not installed')
    existing=tmp_path/'truth.json'
    existing.write_text('caller-owned sentinel')
    result=subprocess.run([str(interpreter),str(REPO/'scripts/smoke_generate.py'),'--output',str(tmp_path)],capture_output=True,text=True)
    assert result.returncode != 0
    assert 'already exists' in result.stderr
    assert existing.read_text()=='caller-owned sentinel'


def test_reference_hashes_include_full_cds_records(tmp_path):
    import subprocess
    interpreter=REPO/'.tool-envs/generation/bin/python'
    if not interpreter.exists():
        pytest.skip('generation runtime not installed')
    refs=tmp_path/'references';refs.mkdir()
    for name in ('P01857-1.fasta','P01834.fasta','J00228.embl','J00241.embl'):
        (refs/name).write_text(name)
    code="""
import hashlib, importlib.util, sys
from pathlib import Path
spec=importlib.util.spec_from_file_location('generation',sys.argv[1]);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
result=mod.reference_hashes(Path(sys.argv[2]))
assert len(result)==4
for name in ('J00228.embl','J00241.embl'):
    p=Path(sys.argv[2])/name
    assert result[str(p)]==hashlib.sha256(p.read_bytes()).hexdigest()
"""
    subprocess.run([str(interpreter),'-c',code,str(REPO/'scripts/smoke_generate.py'),str(refs)],check=True)

@pytest.mark.parametrize('pigz_mode', [False, True])
@pytest.mark.parametrize('read_seed', [2468, 2469])
def test_read_compression_preserves_bytes_and_default_command(tmp_path, monkeypatch, pigz_mode, read_seed):
    from scripts import smoke_reads as m
    import hashlib
    from types import SimpleNamespace
    (tmp_path/'read_counts.tsv').write_text('chain\t2\n')
    pigz=tmp_path/'pigz';pigz.write_bytes(b'constructed native executable');pigz.chmod(0o755)
    content=b'@chain/1\nACGT\n+\nIIII\n'
    commands=[]
    def native(command, **kw):
        commands.append(command)
        assert kw['timeout']==12. and all(kw['env'][k]=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'))
        if command[1]=='generate':
            for suffix in ('R1','R2'):
                raw=tmp_path/f'amplicon_{suffix}.fastq'
                if '--compress' in command:
                    with gzip.open(str(raw)+'.gz','wb') as f:f.write(content)
                else:raw.write_bytes(content)
        else:
            assert command==[str(pigz.resolve()),'-k','-n','-p','3','--',str(tmp_path/'amplicon_R1.fastq'),str(tmp_path/'amplicon_R2.fastq')]
            for raw in map(Path,command[-2:]):
                with gzip.open(str(raw)+'.gz','wb') as f:f.write(raw.read_bytes())
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(m.subprocess,'run',native)
    monkeypatch.setattr(m.subprocess,'check_output',lambda *a,**kw:'pigz 2.4\n')
    record=m.run(tmp_path,3,'amplicon',12.,pigz=pigz if pigz_mode else None,seed=read_seed)[0]
    assert record['pairs']==1 and record['seed']==read_seed
    assert record['logical_fastq']['R1']['sha256']==hashlib.sha256(content).hexdigest()
    assert commands[0][commands[0].index('--seed')+1]==str(read_seed)
    if pigz_mode:
        assert '--compress' not in commands[0] and len(commands)==2
        assert record['compression']['version']=='pigz 2.4'
        for suffix in ('R1','R2'):
            assert record['compression']['logical_fastq'][suffix]==dict(bytes=len(content),sha256=hashlib.sha256(content).hexdigest())
            assert not (tmp_path/f'amplicon_{suffix}.fastq').exists()
    else:
        assert commands==[[str(Path(m.sys.executable).parent/'iss'),'generate','--genomes',str(tmp_path/'amplicon_templates.fasta'),'--readcount_file',str(tmp_path/'read_counts.tsv'),'--sequence_type','amplicon','--model','MiSeq','--seed',str(read_seed),'--cpus','3','--compress','--output',str(tmp_path/'amplicon')]]
        assert 'compression' not in record

@pytest.mark.parametrize('failure', ['native', 'timeout', 'decoded_mismatch'])
def test_pigz_failure_retains_plaintext_and_refuses_existing_output(tmp_path, monkeypatch, failure):
    from scripts import smoke_reads as m
    import subprocess
    (tmp_path/'read_counts.tsv').write_text('chain\t2\n')
    pigz=tmp_path/'pigz';pigz.write_bytes(b'constructed native executable');pigz.chmod(0o755)
    content=b'@chain/1\nACGT\n+\nIIII\n'
    def native(command, **kw):
        if command[1]=='generate':
            for suffix in ('R1','R2'):(tmp_path/f'amplicon_{suffix}.fastq').write_bytes(content)
        else:
            if failure in ('native','timeout'):
                with gzip.open(str(command[-2])+'.gz','wb') as f:f.write(content)
                if failure=='timeout':raise subprocess.TimeoutExpired(command,kw['timeout'])
                raise subprocess.CalledProcessError(1,command)
            for raw in map(Path,command[-2:]):
                with gzip.open(str(raw)+'.gz','wb') as f:f.write(content.replace(b'ACGT',b'TGCA'))
    monkeypatch.setattr(m.subprocess,'run',native)
    monkeypatch.setattr(m.subprocess,'check_output',lambda *a,**kw:'pigz 2.4\n')
    with pytest.raises((subprocess.CalledProcessError,subprocess.TimeoutExpired,ValueError)):
        m.run(tmp_path,3,'amplicon',12.,pigz=pigz)
    for suffix in ('R1','R2'):assert (tmp_path/f'amplicon_{suffix}.fastq').read_bytes()==content
    assert not (tmp_path/'read_run_manifest.json').exists()
    with pytest.raises(FileExistsError):m.run(tmp_path,3,'amplicon',12.,pigz=pigz)


def test_read_manifest_blocks_completed_stage_without_retained_mates(tmp_path, monkeypatch):
    from scripts import smoke_reads as m
    manifest=tmp_path/'read_run_manifest.json'
    manifest.write_text('completed caller-owned receipt')
    def unexpected(*args, **kwargs):raise AssertionError('generation must not start')
    monkeypatch.setattr(m.subprocess,'run',unexpected)
    with pytest.raises(FileExistsError,match='manifest already exists'):
        m.run(tmp_path,3,'amplicon')
    assert manifest.read_text()=='completed caller-owned receipt'


@pytest.mark.parametrize('platform,factor', [('darwin',1),('linux',1024)])
def test_read_child_rss_manifest_uses_bytes(tmp_path, monkeypatch, platform, factor):
    from scripts import smoke_reads as m
    from types import SimpleNamespace
    (tmp_path/'read_counts.tsv').write_text('chain\t2\n')
    def native(command, **kwargs):
        for suffix in ('R1','R2'):
            with gzip.open(tmp_path/f'amplicon_{suffix}.fastq.gz','wb') as f:
                f.write(b'@chain\nACGT\n+\nIIII\n')
    monkeypatch.setattr(m.subprocess,'run',native)
    monkeypatch.setattr(m.sys,'platform',platform)
    monkeypatch.setattr(m.resource,'getrusage',lambda _:SimpleNamespace(ru_maxrss=123))
    assert m.run(tmp_path,3,'amplicon')[0]['peak_child_rss_bytes']==123*factor
