"""Check the small all-locus substrate using installed published generation tools."""
import json
from pathlib import Path
import shutil
import subprocess
import pytest

REPO = Path(__file__).parents[1]


def test_all_loci_paired_internal_branch_and_spliced_constants(tmp_path):
    interpreter = REPO / '.tool-envs/generation/bin/python'
    sources = {name: REPO/'tmp/tool-smoke-complete/references'/name for name in
               ('P01857-1.fasta', 'P01834.fasta', 'J00228.embl', 'J00241.embl')}
    sources.update({name: REPO/'tmp/lambda-reference-audit'/name for name in
                    ('P0DOY2.fasta', 'J00253.embl')})
    if not interpreter.exists() or not all(p.exists() for p in sources.values()):
        pytest.skip('installed generation runtime and verified raw reference records required')
    refs = tmp_path/'references'
    refs.mkdir()
    for name, source in sources.items():
        shutil.copyfile(source, refs/name)
    result = subprocess.run([str(interpreter), str(REPO/'scripts/smoke_generate.py'),
                             '--output', str(tmp_path), '--all-loci'],
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    chains = json.loads((tmp_path/'truth.json').read_text())['chains']
    assert len(chains) == 12
    assert [sum(c['locus'] == locus for c in chains) for locus in ('IGH','IGK','IGL')] == [6,3,3]
    for family in ('family0','family1'):
        members = [c for c in chains if c['family_id'] == family]
        assert len({c['antibody_id'] for c in members}) == 3
        for locus in {c['locus'] for c in members}:
            rows = [c for c in members if c['locus'] == locus]
            assert rows[0]['parent_id'] is None
            assert rows[1]['parent_id'] == rows[0]['chain_id']
            assert rows[2]['parent_id'] == rows[1]['chain_id']
            inherited = [i for i, (a,b) in enumerate(zip(rows[0]['variable_nt'], rows[1]['variable_nt'])) if a != b]
            assert inherited and any(rows[2]['variable_nt'][i] == rows[1]['variable_nt'][i] for i in inherited)
    assert all(c['airr_truth']['j_call'].startswith('IGLJ2/3-') and
               c['constant_nt_gene'] == 'IGLC2' for c in chains if c['locus'] == 'IGL')


def test_read_counts_use_individual_reads_and_reject_odd_counts(tmp_path):
    import importlib.util
    spec = importlib.util.spec_from_file_location('smoke_reads', REPO/'scripts/smoke_reads.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    path = tmp_path/'counts.tsv'
    path.write_text('a\t2\nb\t8\nc\t0\n')
    assert module.expected_pair_count(path) == 5
    path.write_text('a\t3\n')
    with pytest.raises(ValueError):
        module.expected_pair_count(path)


def test_read_generation_selects_only_requested_geometry(tmp_path, monkeypatch):
    import gzip
    import json
    import importlib.util
    spec = importlib.util.spec_from_file_location('selected_reads', REPO/'scripts/smoke_reads.py')
    reads = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reads)
    assert reads.selected_modes() == ('amplicon','rnaseq')
    with pytest.raises(ValueError):
        reads.selected_modes('unknown')
    for geometry in ('amplicon','rnaseq'):
        root=tmp_path/geometry;root.mkdir()
        (root/'read_counts.tsv').write_text('c0\t2\n')
        commands=[]
        def native(command, **kwargs):
            assert kwargs['timeout']==900.
            commands.append(command)
            for suffix in ('R1','R2'):
                with gzip.open(root/f'{geometry}_{suffix}.fastq.gz','wt') as handle:
                    handle.write('@c0\nACG\n+\nIII\n')
        monkeypatch.setattr(reads.subprocess,'run',native)
        records=reads.run(root,3,geometry,stage_timeout=900.)
        assert len(commands)==len(records)==1
        assert records[0]['scenario']==geometry and records[0]['pairs']==1
        assert commands[0][commands[0].index('--sequence_type')+1] == ('amplicon' if geometry=='amplicon' else 'metagenomics')
        assert {r['scenario'] for r in json.loads((root/'read_run_manifest.json').read_text())}=={geometry}
