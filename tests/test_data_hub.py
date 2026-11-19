import hashlib

import pytest

from vdjaqp import data_hub


def _manifest(tmp_path, rows):
    lines = ['path\tbytes\tsha256\tgroup']
    for path, data, group in rows:
        f = tmp_path / 'copy' / path
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_bytes(data)
        lines.append(f'{path}\t{len(data)}\t{hashlib.sha256(data).hexdigest()}\t{group}')
    m = tmp_path / 'MANIFEST.tsv'
    m.write_text('\n'.join(lines) + '\n')
    return m


def test_select_defaults_to_calibration(tmp_path):
    m = _manifest(tmp_path, [('calibration/model.json', b'{}', 'calibration'),
                             ('results/a.parquet', b'x' * 10, 'results')])
    rows = data_hub.read_manifest(m)
    assert [r['path'] for r in data_hub.select(rows, data_hub.DEFAULT_GROUPS)] == ['calibration/model.json']
    assert len(data_hub.select(rows, ['all'])) == 2


def test_verify_detects_corruption(tmp_path):
    m = _manifest(tmp_path, [('calibration/model.json', b'{"a":1}', 'calibration')])
    rows = data_hub.read_manifest(m)
    assert data_hub.verify(rows, tmp_path / 'copy') == []
    (tmp_path / 'copy/calibration/model.json').write_bytes(b'{"a":2}')
    assert data_hub.verify(rows, tmp_path / 'copy') == ['sha256 calibration/model.json']


def test_unknown_group_rejected(tmp_path):
    m = tmp_path / 'MANIFEST.tsv'
    m.write_text('path\tbytes\tsha256\tgroup\nx\t1\t00\tprivate\n')
    with pytest.raises(ValueError):
        data_hub.read_manifest(m)


def test_dry_run(tmp_path, capsys):
    m = _manifest(tmp_path, [('calibration/model.json', b'{}', 'calibration')])
    assert data_hub.main(['--dry-run', '--manifest', str(m)]) == 0
    assert '1 files' in capsys.readouterr().out
