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


def test_bootstrap_and_library_share_portable_root(tmp_path, monkeypatch):
    import importlib
    from pathlib import Path
    from vdjaqp import paths
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv('VDJAQP_DATA', raising=False)
    destinations = []
    monkeypatch.setattr(data_hub, 'fetch', lambda groups, dest, revision: destinations.append(dest))
    try:
        for configured in (None, '~/vdjaqp-example'):
            if configured is not None: monkeypatch.setenv('VDJAQP_DATA', configured)
            importlib.reload(paths)
            assert data_hub.main([]) == 0
            expected = Path(configured or 'data').expanduser()
            assert destinations[-1] == paths.DATA == expected
            assert paths.GUT == expected/'external'/'gut'
    finally:
        monkeypatch.delenv('VDJAQP_DATA', raising=False)
        importlib.reload(paths)

@pytest.mark.parametrize('content', [
    '', 'path\tgroup\n',
    'path\tbytes\tsha256\tgroup\n../escape\t0\t' + '0'*64 + '\tcalibration\n',
    'path\tbytes\tsha256\tgroup\n/escape\t0\t' + '0'*64 + '\tcalibration\n',
    'path\tbytes\tsha256\tgroup\nx\t-1\t' + '0'*64 + '\tcalibration\n',
    'path\tbytes\tsha256\tgroup\nx\t0\tbad\tcalibration\n',
    'path\tbytes\tsha256\tgroup\nx\t0\t' + '0'*64 + '\tcalibration\nx\t1\t' + '1'*64 + '\tcalibration\n',
])
def test_manifest_rejects_invalid_metadata(tmp_path, content):
    manifest = tmp_path/'MANIFEST.tsv'
    manifest.write_text(content)
    with pytest.raises(ValueError):
        data_hub.read_manifest(manifest)


def _fake_hub(monkeypatch, snapshot):
    import sys
    import types
    revisions = []
    def download(*args, revision=None, **kwargs):
        revisions.append(revision)
        return str(snapshot/'MANIFEST.tsv')
    def snapshot_download(*args, revision=None, **kwargs):
        revisions.append(revision)
        return str(snapshot)
    monkeypatch.setitem(sys.modules, 'huggingface_hub', types.SimpleNamespace(
        HfApi=lambda: types.SimpleNamespace(dataset_info=lambda *a, **k: types.SimpleNamespace(sha='a'*40)),
        hf_hub_download=download, snapshot_download=snapshot_download))
    return revisions


def _snapshot(tmp_path, payload=b'new'):
    m = _manifest(tmp_path, [('calibration/model.json', payload, 'calibration')])
    snapshot = tmp_path/'copy'
    (snapshot/'MANIFEST.tsv').write_bytes(m.read_bytes())
    return snapshot


def test_fetch_pins_one_revision_and_rejects_destination_escape(tmp_path, monkeypatch):
    snapshot = _snapshot(tmp_path)
    revisions = _fake_hub(monkeypatch, snapshot)
    outside = tmp_path/'outside'
    outside.mkdir()
    (outside/'model.json').write_bytes(b'old')
    dest = tmp_path/'dest'
    dest.mkdir()
    (dest/'calibration').symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match='destination'):
        data_hub.fetch(['calibration'], dest, revision='main')
    assert (outside/'model.json').read_bytes() == b'old'
    assert revisions == ['a'*40, 'a'*40]


def test_interrupted_copy_preserves_existing_file(tmp_path, monkeypatch):
    snapshot = _snapshot(tmp_path)
    _fake_hub(monkeypatch, snapshot)
    dest = tmp_path/'dest'
    target = dest/'calibration/model.json'
    target.parent.mkdir(parents=True)
    target.write_bytes(b'accepted')
    def interrupted(source, target):
        from pathlib import Path
        Path(target).write_bytes(b'partial')
        raise OSError('interrupted')
    monkeypatch.setattr(data_hub.shutil, 'copyfile', interrupted)
    with pytest.raises(OSError, match='interrupted'):
        data_hub.fetch(['calibration'], dest)
    assert target.read_bytes() == b'accepted'
    assert list(target.parent.iterdir()) == [target]


def test_bad_staged_copy_never_replaces_accepted_bytes(tmp_path, monkeypatch):
    snapshot = _snapshot(tmp_path)
    _fake_hub(monkeypatch, snapshot)
    dest = tmp_path/'dest'
    target = dest/'calibration/model.json'
    target.parent.mkdir(parents=True)
    target.write_bytes(b'accepted')
    monkeypatch.setattr(data_hub.shutil, 'copyfile', lambda source, target: target.write_bytes(b'bad'))
    with pytest.raises(ValueError, match='staged'):
        data_hub.fetch(['calibration'], dest)
    assert target.read_bytes() == b'accepted'


def test_local_lock_bounds_wait_and_recovers_after_process_exit(tmp_path):
    import subprocess
    import sys
    import time
    ready = tmp_path/'ready'
    code = '''
import sys, time
from pathlib import Path
from vdjaqp.data_hub import _destination_lock
root = Path(sys.argv[1])
with _destination_lock(root):
    (root/'ready').touch()
    time.sleep(20)
'''
    owner = subprocess.Popen([sys.executable, '-c', code, str(tmp_path)])
    try:
        deadline = time.monotonic() + 10
        while not ready.exists() and time.monotonic() < deadline:
            assert owner.poll() is None
            time.sleep(0.02)
        assert ready.exists()
        with pytest.raises(TimeoutError):
            with data_hub._destination_lock(tmp_path, timeout=0.05):
                pytest.fail('entered another process lock')
        owner.kill()
        owner.wait(timeout=5)
        with data_hub._destination_lock(tmp_path, timeout=1):
            pass
    finally:
        if owner.poll() is None:
            owner.kill()
        owner.wait(timeout=5)


def test_concurrent_fetches_serialize_complete_writers(tmp_path):
    import subprocess
    import sys
    import time
    first = _snapshot(tmp_path/'first', b'first')
    second = _snapshot(tmp_path/'second', b'second')
    dest, ready, release = (tmp_path/name for name in ('dest', 'ready', 'release'))
    code = '''
import sys, types, time
from pathlib import Path
from vdjaqp import data_hub
snapshot, dest, ready, release = map(Path, sys.argv[1:5])
first = sys.argv[5] == 'first'
sys.modules['huggingface_hub'] = types.SimpleNamespace(
    HfApi=lambda: types.SimpleNamespace(dataset_info=lambda *a, **k: types.SimpleNamespace(sha='a'*40)),
    hf_hub_download=lambda *a, **k: str(snapshot/'MANIFEST.tsv'),
    snapshot_download=lambda *a, **k: str(snapshot))
copy = data_hub.shutil.copyfile
if first:
    def delayed(source, target):
        ready.touch()
        deadline = time.monotonic()+10
        while not release.exists():
            if time.monotonic()>deadline: raise TimeoutError('test barrier')
            time.sleep(0.02)
        return copy(source, target)
    data_hub.shutil.copyfile = delayed
data_hub.fetch(['calibration'], dest)
'''
    processes = []
    try:
        processes.append(subprocess.Popen([sys.executable, '-c', code, str(first), str(dest), str(ready), str(release), 'first'], stdout=subprocess.DEVNULL))
        deadline = time.monotonic()+10
        while not ready.exists() and time.monotonic()<deadline:
            assert processes[0].poll() is None
            time.sleep(0.02)
        assert ready.exists()
        processes.append(subprocess.Popen([sys.executable, '-c', code, str(second), str(dest), str(ready), str(release), 'second'], stdout=subprocess.DEVNULL))
        with pytest.raises(subprocess.TimeoutExpired):
            processes[1].wait(timeout=0.2)
        release.touch()
        for process in processes:
            assert process.wait(timeout=10) == 0
        assert (dest/'calibration/model.json').read_bytes() == b'second'
        assert data_hub.verify(data_hub.read_manifest(dest/'MANIFEST.igprotseq.tsv'), dest) == []
    finally:
        release.touch()
        for process in processes:
            if process.poll() is None:
                process.kill()
            process.wait(timeout=5)


@pytest.mark.parametrize('name', ['a\\b', 'a/../b', './a', 'a//b', '.vdjaqp-data.lock', 'MANIFEST.igprotseq.tsv'])
def test_manifest_path_cannot_alias_or_replace_metadata(name):
    with pytest.raises(ValueError):
        data_hub._relative_path(name)


@pytest.mark.parametrize('name', ['C:/escape', 'C:escape', '//server/share/file', '\\\\server\\share\\file'])
def test_manifest_rejects_windows_path_roots(name):
    with pytest.raises(ValueError):
        data_hub._relative_path(name)


def test_lock_preserves_inode_and_rejects_symlinks(tmp_path):
    lock = tmp_path/'.vdjaqp-data.lock'
    with data_hub._destination_lock(tmp_path):
        identity = lock.stat().st_ino
    assert lock.stat().st_ino == identity
    with data_hub._destination_lock(tmp_path):
        assert lock.stat().st_ino == identity
    lock.unlink()
    outside = tmp_path/'outside'
    outside.write_bytes(b'accepted')
    lock.symlink_to(outside)
    with pytest.raises(ValueError, match='symlink'):
        with data_hub._destination_lock(tmp_path):
            pytest.fail('locked symlink')
    assert outside.read_bytes() == b'accepted'


@pytest.mark.parametrize('alias_kind', ['directory_symlink', 'hardlink', 'case_alias'])
def test_fetch_cannot_replace_held_lock_through_alias(tmp_path, monkeypatch, alias_kind):
    dest = tmp_path/'dest'
    dest.mkdir()
    lock = dest/'.vdjaqp-data.lock'
    lock.touch()
    identity = lock.stat().st_ino
    if alias_kind == 'directory_symlink':
        (dest/'alias').symlink_to(dest, target_is_directory=True)
        malicious = 'alias/.vdjaqp-data.lock'
    elif alias_kind == 'hardlink':
        (dest/'lock-alias').hardlink_to(lock)
        malicious = 'lock-alias'
    else:
        malicious = '.VDJAQP-DATA.LOCK'
        if not (dest/malicious).exists():
            pytest.skip('case-sensitive filesystem has no case alias')
    manifest = _manifest(tmp_path, [('calibration/model.json', b'new', 'calibration'),
                                     (malicious, b'replacement', 'calibration')])
    snapshot = tmp_path/'copy'
    (snapshot/'MANIFEST.tsv').write_bytes(manifest.read_bytes())
    _fake_hub(monkeypatch, snapshot)
    with pytest.raises(ValueError, match='bootstrap lock'):
        data_hub.fetch(['calibration'], dest)
    assert lock.stat().st_ino == identity
    assert not (dest/'calibration/model.json').exists()
    with data_hub._destination_lock(dest, timeout=1):
        pass
