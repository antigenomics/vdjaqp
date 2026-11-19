"""The documented conda install must verify the existing artifact pins."""
import json
from pathlib import Path


def test_macos_explicit_lock_uses_manifest_hashes():
    root=Path(__file__).parents[1]/'envs/preprocess'
    packages=json.loads((root/'runtime.json').read_text())['packages']
    entries=[line for line in (root/'osx-64-explicit.txt').read_text().splitlines()
        if line.startswith('https://')]
    assert set(entries)=={p['url']+'#'+p['md5'] for p in packages}
