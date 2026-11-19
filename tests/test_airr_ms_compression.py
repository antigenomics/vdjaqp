"""Lossless final-artifact packaging and pre-mutation refusal boundaries."""
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import pytest

spec = importlib.util.spec_from_file_location('pack', Path(__file__).parents[1]/'scripts/compress_airr_ms_artifacts.py')
pack = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pack)


def test_lossless_completed_case_packaging_and_refusals(tmp_path, monkeypatch):
    monkeypatch.setattr(pack, 'MIN_BYTES', 10)
    payload = b'private final table\n' * 20
    def setup(name):
        root = tmp_path/name
        root.mkdir()
        (root/'final.tsv').write_bytes(payload)
        (root/'keep.parquet').write_bytes(payload)
        (root/'small.tsv').write_bytes(b'small')
        report = dict(spectra_cleanup={'run/spectra.mzML': {'status': 'deleted'}},
            artifact_sha256={f: pack.digest(root/f) for f in ('final.tsv','keep.parquet','small.tsv')})
        (root/'airr_ms_report.json').write_text(json.dumps(report))
        return root, report
    root, report = setup('good')
    before = (root/'airr_ms_report.json').read_bytes()
    retention = pack.compress_case(root)
    metadata = json.loads(retention.read_text())
    assert metadata['status'] == 'complete'
    entry, = metadata['artifacts']
    assert entry['status'] == 'compressed_original_removed'
    assert gzip.decompress((root/'final.tsv.gz').read_bytes()) == payload
    assert entry['original_sha256'] == hashlib.sha256(payload).hexdigest()
    assert entry['gzip_sha256'] == pack.digest(root/'final.tsv.gz')
    assert not (root/'final.tsv').exists()
    assert (root/'keep.parquet').read_bytes() == payload
    assert (root/'small.tsv').read_bytes() == b'small'
    assert (root/'airr_ms_report.json').read_bytes() == before
    with pytest.raises(FileExistsError):
        pack.compress_case(root)
    for mode in ('tampered','pending','escape','existing_gzip'):
        root, report = setup(mode)
        if mode == 'tampered': report['artifact_sha256']['final.tsv'] = 'bad'
        if mode == 'pending': report['spectra_cleanup']['run/spectra.mzML']['status'] = 'pending'
        if mode == 'escape': report['artifact_sha256']['../outside.tsv'] = 'bad'
        if mode == 'existing_gzip': (root/'final.tsv.gz').write_bytes(b'existing')
        (root/'airr_ms_report.json').write_text(json.dumps(report))
        with pytest.raises((ValueError, FileExistsError)):
            pack.compress_case(root)
        assert (root/'final.tsv').read_bytes() == payload
        assert not (root/'retention.json').exists()
