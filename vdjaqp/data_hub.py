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
import shutil
import sys
from pathlib import Path

REPO_ID = 'isalgo/igprotseq'
REPO_TYPE = 'dataset'
GROUPS = ('calibration', 'donor', 'simulations', 'external', 'results')
DEFAULT_GROUPS = ('calibration',)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read_manifest(path: Path) -> list[dict]:
    with path.open(newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    for r in rows:
        if not {'path', 'bytes', 'sha256', 'group'} <= r.keys():
            raise ValueError(f'manifest row lacks required columns: {r}')
        if r['group'] not in GROUPS:
            raise ValueError(f"unknown group {r['group']!r} for {r['path']}")
    return rows


def select(rows: list[dict], groups) -> list[dict]:
    groups = set(GROUPS if 'all' in groups else groups)
    return [r for r in rows if r['group'] in groups]


def verify(rows: list[dict], root: Path, strip: str = '') -> list[str]:
    """Return a list of problems (missing, size or sha256 mismatch) under ``root``."""
    problems = []
    for r in rows:
        rel = r['path'][len(strip):] if strip and r['path'].startswith(strip) else r['path']
        p = root / rel
        if not p.is_file():
            problems.append(f'missing {rel}')
        elif p.stat().st_size != int(r['bytes']):
            problems.append(f'size {rel}')
        elif sha256(p) != r['sha256']:
            problems.append(f'sha256 {rel}')
    return problems


def fetch(groups, dest: Path, revision: str | None = None, cache: Path | None = None) -> Path:
    from huggingface_hub import hf_hub_download, snapshot_download
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
    dest.mkdir(parents=True, exist_ok=True)
    for r in rows:
        # The results group mirrors the paper repository layout: results/... -> DEST/results/...
        target = dest / r['path']
        target.parent.mkdir(parents=True, exist_ok=True)
        src = snapshot / r['path']
        if target.exists() and target.stat().st_size == int(r['bytes']) and sha256(target) == r['sha256']:
            continue
        shutil.copyfile(src, target)
    shutil.copyfile(manifest, dest / 'MANIFEST.igprotseq.tsv')
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
