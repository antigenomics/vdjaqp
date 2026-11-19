"""Reference integrity checks only; no repertoire generation or network access."""
import json
from pathlib import Path
import shutil
import subprocess

import pytest


REPO = Path(__file__).resolve().parents[1]


def test_constant_candidates_and_tampered_input(tmp_path):
    python = REPO / '.tool-envs/generation/bin/python'
    raw = REPO / 'tmp/constant-reference-audit'
    manifest = json.loads((REPO / 'envs/generation/constant_candidates.json').read_text())
    if not python.exists() or not all((raw / n).exists() for n in manifest['downloaded_records']):
        pytest.skip('generation runtime and pinned raw constant reference inputs required')
    for name in manifest['downloaded_records']:
        shutil.copyfile(raw / name, tmp_path / name)
    command = [str(python), str(REPO / 'scripts/validate_constant_candidates.py'),
               '--reference-dir', str(tmp_path)]
    checked = subprocess.run(command, capture_output=True, text=True, timeout=15)
    assert checked.returncode == 0, checked.stderr
    rows = json.loads(checked.stdout)
    assert [(r['gene'], r['physical_constant_aa']) for r in rows] == [
        ('IGHM', 453), ('IGHG3', 362), ('IGHA1', 353), ('IGHE', 428)]
    # The pinned runtime owns Biopython; pytest itself needs no extra dependency.
    api_check = '''
import hashlib
import json
import sys
from pathlib import Path
from Bio.Seq import Seq
from scripts.validate_constant_candidates import load_validated_constants, join_constant
root = Path(sys.argv[1])
before = {p.name: p.read_bytes() for p in root.iterdir()}
records = load_validated_constants(root)
manifest = json.loads(Path('envs/generation/constant_candidates.json').read_text())
assert list(records) == [r['gene'] for r in manifest['candidates']]
for expected in manifest['candidates']:
    r = records[expected['gene']]
    assert r['cds_accession'] == expected['nucleotide_accession']
    assert r['protein_accession'] == expected['protein_id']
    assert hashlib.sha256(r['raw_cds_nt'].encode()).hexdigest() == r['cds_sha256'] == expected['cds_sha256']
    assert hashlib.sha256(r['aa'].encode()).hexdigest() == expected['physical_constant_sha256']
    assert r['stop_codon'] == r['raw_cds_nt'][-3:] in {'TAA', 'TAG', 'TGA'}
    assert r['nt'] == (r['raw_cds_nt'][1:] if r['gene'] in {'IGHM', 'IGHG3'} else r['raw_cds_nt'])
    joined = join_constant('ATGG', 'M', r)
    assert joined == 'ATGG' + r['nt']
    assert str(Seq(joined).translate()) == 'M' + r['aa'] + '*'
    for bad in ('ATGA', 'ATGC', 'ATGT', 'ATG'):
        try:
            join_constant(bad, 'M', r)
        except ValueError:
            pass
        else:
            raise AssertionError('incompatible native tail accepted')
    try:
        join_constant('ATGG', 'L', r)
    except ValueError:
        pass
    else:
        raise AssertionError('incorrect variable translation accepted')
assert before == {p.name: p.read_bytes() for p in root.iterdir()}
'''
    checked_api = subprocess.run([str(python), '-c', api_check, str(tmp_path)],
                                 cwd=REPO, capture_output=True, text=True, timeout=15)
    assert checked_api.returncode == 0, checked_api.stderr
    (tmp_path / 'P01871-1.fasta').write_text('>changed\nG\n')
    failed = subprocess.run(command, capture_output=True, text=True, timeout=15)
    assert failed.returncode != 0
    assert 'raw SHA256 mismatch' in failed.stderr
