"""The integration driver must not overwrite caller-selected directories."""
from pathlib import Path
import subprocess
import sys


def test_existing_output_is_refused_before_tools_or_network(tmp_path):
    marker = tmp_path / 'keep.txt'
    marker.write_text('preserve')
    driver = Path(__file__).resolve().parents[1] / 'scripts/smoke_all.py'
    result = subprocess.run([sys.executable, str(driver), '--output', str(tmp_path)],
                            capture_output=True, text=True, timeout=10)
    assert result.returncode != 0 and 'already exists' in result.stderr
    assert marker.read_text() == 'preserve'
def test_raw_reference_input_is_checksum_verified_without_network(tmp_path):
    import hashlib
    import importlib.util
    import json
    from pathlib import Path
    import pytest
    script = Path(__file__).parents[1]/'scripts/smoke_all.py'
    spec = importlib.util.spec_from_file_location('raw_reference_driver', script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    repo = tmp_path/'repo'
    (repo/'envs/generation').mkdir(parents=True)
    raw = tmp_path/'raw'
    raw.mkdir()
    content = b'>public_reference\nACGT\n'
    (raw/'record.fasta').write_bytes(content)
    manifest = {'references':{'record.fasta':{'url':'invalid://must-not-fetch',
                 'sha256':hashlib.sha256(content).hexdigest()}}}
    (repo/'envs/generation/references.json').write_text(json.dumps(manifest))
    module.prepare_references(repo,tmp_path/'run',source_dir=raw)
    assert (tmp_path/'run/references/record.fasta').read_bytes() == content
    (raw/'record.fasta').write_bytes(b'corrupted')
    with pytest.raises(ValueError,match='hash mismatch'):
        module.prepare_references(repo,tmp_path/'bad-run',source_dir=raw)
