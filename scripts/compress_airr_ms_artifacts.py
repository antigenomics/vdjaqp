"""Losslessly package listed final AIRR-MS text artifacts after completed cleanup.

No scoring or reconstruction runs; original result reports remain unchanged.
Interrupted packaging retains its separate manifest for explicit recovery.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile

MIN_BYTES = 1_000_000
SUFFIXES = {'.tsv', '.idXML', '.featureXML', '.json', '.fasta'}


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def contained(root, name):
    relative = Path(name)
    if relative.is_absolute() or '..' in relative.parts or not relative.parts:
        raise ValueError('artifact path must remain relative to the completed case')
    path = root / relative
    if not path.resolve().is_relative_to(root) or any(p.is_symlink() for p in [path, *path.parents] if p != root and p.is_relative_to(root)):
        raise ValueError('artifact path escapes case or traverses a symlink')
    return path


def compress_case(case):
    root = Path(case).resolve(strict=True)
    report_path = contained(root, 'airr_ms_report.json')
    retention = root / 'retention.json'
    if retention.exists() or retention.is_symlink():
        raise FileExistsError('retention manifest already exists; inspect prior packaging')
    report = json.loads(report_path.read_text())
    cleanup = report.get('spectra_cleanup')
    if not cleanup or any(item.get('status') != 'deleted' or contained(root, name).exists()
                          for name, item in cleanup.items()):
        raise ValueError('case must have completed native spectra cleanup')
    selected = []
    for name, expected in report['artifact_sha256'].items():
        raw = contained(root, name)
        if raw.suffix not in SUFFIXES:
            continue
        if not raw.is_file():
            raise ValueError('listed text artifact is missing')
        size = raw.stat().st_size
        if size < MIN_BYTES:
            continue
        compressed = raw.with_name(raw.name + '.gz')
        if compressed.exists() or compressed.is_symlink():
            raise FileExistsError('compressed artifact already exists')
        if digest(raw) != expected:
            raise ValueError('raw artifact SHA256 differs from final report')
        selected.append((name, raw, compressed, expected, size))
    manifest = dict(report_sha256=digest(report_path), status='in_progress', artifacts=[])
    # Exclusive creation reserves this packaging operation; never resume silently.
    with retention.open('x') as stream:
        json.dump(manifest, stream, indent=2)
    def save():
        with tempfile.NamedTemporaryFile(mode='w', dir=root, prefix='.retention-', delete=False) as stream:
            temp = Path(stream.name)
            json.dump(manifest, stream, indent=2)
            stream.write('\n')
        temp.replace(retention)
    for name, raw, compressed, expected, size in selected:
        with tempfile.NamedTemporaryFile(dir=raw.parent, prefix='.'+raw.name+'-', delete=False) as stream:
            temporary = Path(stream.name)
        try:
            with raw.open('rb') as source, temporary.open('wb') as target:
                with gzip.GzipFile(filename='', fileobj=target, mode='wb', mtime=0) as zipped:
                    shutil.copyfileobj(source, zipped, length=1024*1024)
            with gzip.open(temporary, 'rb') as restored:
                if hashlib.file_digest(restored, 'sha256').hexdigest() != expected:
                    raise ValueError('gzip roundtrip SHA256 mismatch')
            entry = dict(raw_path=name, original_sha256=expected, raw_bytes=size,
                         gzip_path=str(compressed.relative_to(root)), gzip_sha256=digest(temporary),
                         gzip_bytes=temporary.stat().st_size, status='verified_not_published')
            manifest['artifacts'].append(entry)
            save()
            # link is an atomic, no-overwrite publication on the sibling filesystem.
            os.link(temporary, compressed)
            entry['status'] = 'published_original_retained'
            save()
            raw.unlink()
            entry['status'] = 'compressed_original_removed'
            save()
        finally:
            temporary.unlink(missing_ok=True)
    manifest['status'] = 'complete'
    save()
    return retention


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', type=Path)
    print(compress_case(parser.parse_args().case))
