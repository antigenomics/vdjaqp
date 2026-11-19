"""Fetch VDJ-AQP data from the Hugging Face dataset ``isalgo/igprotseq``.

The dataset is the single data source for the pipeline and the manuscript. Its
``MANIFEST.tsv`` lists every file with ``path``, ``bytes``, ``sha256`` and ``group``.
By default only the small ``calibration`` group (frozen error model and declared
parameters) is downloaded; larger groups are fetched on request:

    vdjaqp-data                          # calibration only
    vdjaqp-data --group donor simulations
    vdjaqp-data --group results --dest ../2026-vdjaqp-ms   # paper result tables
    vdjaqp-data --group all
    vdjaqp-data --verify DIR             # re-check sha256 of an existing copy

Access follows the dataset settings; use ``hf auth login`` when required.
Every downloaded file is verified against the manifest sha256.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
from contextlib import contextmanager
from tempfile import NamedTemporaryFile
import shutil
import sys
from pathlib import Path, PurePosixPath, PureWindowsPath

REPO_ID = 'isalgo/igprotseq'
REPO_TYPE = 'dataset'
GROUPS = ('calibration', 'donor', 'simulations', 'external', 'results')
DEFAULT_GROUPS = ('calibration',)


def sha256(path: Path) -> str:
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def _relative_path(value: str) -> str:
    if (not isinstance(value, str) or not value or '\\' in value or
            any(ord(c) < 32 for c in value) or
            PurePosixPath(value).is_absolute() or PureWindowsPath(value).drive or
            any(part in ('', '.', '..') for part in value.split('/'))):
        raise ValueError('manifest path must be a normalized relative POSIX path')
    if value.split('/')[0] in ('.vdjaqp-data.lock', 'MANIFEST.igprotseq.tsv'):
        raise ValueError('manifest path conflicts with bootstrap metadata')
    return value


def read_manifest(path: Path) -> list[dict]:
    with path.open(newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames or []
        if (not {'path', 'bytes', 'sha256', 'group'} <= set(fields) or
                len(fields) != len(set(fields))):
            raise ValueError('manifest lacks unique required columns')
        rows = list(reader)
    seen = set()
    for number, r in enumerate(rows, 2):
        # Do not expose malformed rows: manifests may contain private metadata.
        if (None in r or any(r.get(k) is None for k in fields) or
                r['group'] not in GROUPS or
                not re.fullmatch(r'[0-9]+', r['bytes']) or
                not re.fullmatch(r'[0-9a-fA-F]{64}', r['sha256'])):
            raise ValueError(f'invalid manifest metadata at row {number}')
        name = _relative_path(r['path'])
        if name in seen:
            raise ValueError(f'duplicate manifest path at row {number}')
        seen.add(name)
        r['sha256'] = r['sha256'].lower()
    return rows


def _destination(root: Path, relative: str) -> Path:
    target = root / relative
    if not target.resolve().is_relative_to(root.resolve()):
        raise ValueError('destination path escapes the selected root')
    return target


@contextmanager
def _destination_lock(root: Path, timeout: float = 30):
    """Serialize local writers with a persistent, bounded kernel lock.

    Shared-filesystem semantics need qualification. Preserve the lock pathname
    across release; filelock refuses soft-lock fallback for this mode.
    """
    from filelock import FileLock
    lock = root / '.vdjaqp-data.lock'
    if lock.is_symlink():
        raise ValueError('destination lock must not be a symlink')
    with FileLock(lock, timeout=timeout, mode=0o600,
                  preserve_lock_file=True, fallback_to_soft=False):
        yield


def _publish(source: Path, target: Path, size: int, digest: str) -> None:
    """Validate an owned same-filesystem staging file before atomic publication."""
    temporary = None
    try:
        with NamedTemporaryFile(dir=target.parent, prefix='.vdjaqp-', delete=False) as stream:
            temporary = Path(stream.name)
        shutil.copyfile(source, temporary)
        if temporary.stat().st_size != size or sha256(temporary) != digest:
            raise ValueError('staged dataset file failed verification')
        with temporary.open('rb') as stream:
            os.fsync(stream.fileno())
        os.replace(temporary, target)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def select(rows: list[dict], groups) -> list[dict]:
    groups = set(GROUPS if 'all' in groups else groups)
    return [r for r in rows if r['group'] in groups]


def verify(rows: list[dict], root: Path, strip: str = '') -> list[str]:
    """Return a list of problems (missing, size or sha256 mismatch) under ``root``."""
    problems = []
    for r in rows:
        rel = r['path'][len(strip):] if strip and r['path'].startswith(strip) else r['path']
        _relative_path(rel)
        p = root / rel
        if not p.is_file():
            problems.append(f'missing {rel}')
        elif p.stat().st_size != int(r['bytes']):
            problems.append(f'size {rel}')
        elif sha256(p) != r['sha256']:
            problems.append(f'sha256 {rel}')
    return problems


def fetch(groups, dest: Path, revision: str | None = None, cache: Path | None = None) -> Path:
    from huggingface_hub import HfApi, hf_hub_download, snapshot_download
    # Resolve mutable tags/branches once before either download.
    revision = HfApi().dataset_info(REPO_ID, revision=revision).sha
    manifest = Path(hf_hub_download(REPO_ID, 'MANIFEST.tsv', repo_type=REPO_TYPE, revision=revision,
                                    cache_dir=cache))
    rows = select(read_manifest(manifest), groups)
    if not rows:
        sys.exit(f'no files in groups {groups}')
    snapshot = Path(snapshot_download(REPO_ID, repo_type=REPO_TYPE, revision=revision, cache_dir=cache,
                                      allow_patterns=[r['path'] for r in rows] + ['MANIFEST.tsv', 'README.md']))
    problems = verify(rows, snapshot)
    if problems:
        sys.exit('verification failed:\n  ' + '\n  '.join(problems[:20]))
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    with _destination_lock(dest):
        # Validate all destinations before replacing any existing asset.
        targets = [_destination(dest, r['path']) for r in rows]
        manifest_target = _destination(dest, 'MANIFEST.igprotseq.tsv')
        for r, target in zip(rows, targets):
            target.parent.mkdir(parents=True, exist_ok=True)
            _destination(dest, r['path'])
            if target.is_file() and target.stat().st_size == int(r['bytes']) and sha256(target) == r['sha256']:
                continue
            _publish(snapshot / r['path'], target, int(r['bytes']), r['sha256'])
        _publish(manifest, manifest_target, manifest.stat().st_size, sha256(manifest))
    print(f'{len(rows)} files from {REPO_ID} ({", ".join(sorted(set(r["group"] for r in rows)))}) '
          f'verified into {dest}')
    return dest


def main(argv=None):
    from .paths import data_root
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--group', nargs='+', default=list(DEFAULT_GROUPS), choices=[*GROUPS, 'all'])
    p.add_argument('--dest', type=Path, default=data_root(),
                   help='destination root (default: $VDJAQP_DATA or ./data)')
    p.add_argument('--revision', help='dataset revision (tag or commit) to pin')
    p.add_argument('--verify', type=Path, metavar='DIR', help='verify an existing copy against MANIFEST.tsv')
    p.add_argument('--manifest', type=Path, help='local MANIFEST.tsv for --verify or --dry-run')
    p.add_argument('--dry-run', action='store_true', help='list the selected files without downloading')
    a = p.parse_args(argv)
    if a.verify or a.dry_run:
        manifest = a.manifest or (a.verify / 'MANIFEST.igprotseq.tsv' if a.verify else None)
        if manifest is None:
            p.error('--dry-run needs --manifest')
        rows = select(read_manifest(manifest), a.group)
        if a.dry_run:
            print(f'{len(rows)} files, {sum(int(r["bytes"]) for r in rows) / 1e6:.1f} MB in groups {a.group}')
            return 0
        problems = verify(rows, a.verify)
        print('OK' if not problems else '\n'.join(problems))
        return 1 if problems else 0
    fetch(a.group, a.dest, a.revision)
    return 0


if __name__ == '__main__':
    sys.exit(main())
