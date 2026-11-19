#!/bin/sh
# Reproduce the existing pinned source-overlay recipe; never modify the installed wheel.
set -eu
export IGSEQ_ROOT="$PWD"
export ARDA_STAGE="$IGSEQ_ROOT/.tool-envs/arda-overlay-read-map"
test ! -e "$ARDA_STAGE"
mkdir -p "$ARDA_STAGE"
uv pip install --python "$IGSEQ_ROOT/.venv/bin/python" --no-deps \
  --target "$ARDA_STAGE/baseline" arda-mapper==2.36.0 seqtree==1.1.0
git clone --depth 1 --branch v2.36.0 --no-checkout \
  https://github.com/antigenomics/arda.git "$ARDA_STAGE/source"
git -C "$ARDA_STAGE/source" checkout --detach \
  049553f35114f2578c99355f8afd4c107c42c4d9
test "$(git -C "$ARDA_STAGE/source" rev-parse HEAD)" = \
  049553f35114f2578c99355f8afd4c107c42c4d9
"$IGSEQ_ROOT/.venv/bin/python" - <<'PY'
import hashlib, json, os, subprocess
from pathlib import Path
root=Path(os.environ['IGSEQ_ROOT']);stage=Path(os.environ['ARDA_STAGE'])
recipe=root/'envs/arda';source=stage/'source'
m=json.loads((recipe/'source-manifest.json').read_text())
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
patch=recipe/m['patch']['file']
assert patch.stat().st_size==m['patch']['bytes'] and sha(patch)==m['patch']['sha256']
assert sha(recipe/m['upstream_license']['file'])==m['upstream_license']['sha256']
assert sha(source/'LICENSE')==m['upstream_license']['sha256']
for name,pin in m['modules'].items():
    assert sha(source/'src/arda'/name)==pin['baseline_sha256'],name
subprocess.run(['git','apply','--check',str(patch)],cwd=source,check=True)
subprocess.run(['git','apply',str(patch)],cwd=source,check=True)
extra=recipe/m['read_membership_patch']['file']
assert extra.stat().st_size==m['read_membership_patch']['bytes'] and sha(extra)==m['read_membership_patch']['sha256']
subprocess.run(['git','apply','--check',str(extra)],cwd=source,check=True)
subprocess.run(['git','apply',str(extra)],cwd=source,check=True)
for name,pin in m['modules'].items():
    assert sha(source/'src/arda'/name)==pin['patched_sha256'],name
PY
"$IGSEQ_ROOT/.venv/bin/python" - <<'PY'
import hashlib, json, os, shutil
from importlib.metadata import Distribution
from pathlib import Path
root=Path(os.environ['IGSEQ_ROOT']);stage=Path(os.environ['ARDA_STAGE'])
m=json.loads((root/'envs/arda/source-manifest.json').read_text())
baseline=stage/'baseline'
metadata,=baseline.glob('arda_mapper-*.dist-info')
assert Distribution.at(metadata).version=='2.36.0'
installed=baseline/'arda'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
for name,pin in m['modules'].items():
    assert sha(installed/name)==pin['baseline_sha256'],name
before={str(p.relative_to(installed)):sha(p) for p in installed.rglob('*')
        if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ('.pyc','.pyo')}
overlay=stage/'overlay';(overlay/'arda').mkdir(parents=True,exist_ok=False)
shutil.copytree(metadata, overlay/metadata.name)
shutil.copytree(baseline/'seqtree', overlay/'seqtree')
seqtree_metadata,=baseline.glob('seqtree-*.dist-info')
shutil.copytree(seqtree_metadata, overlay/seqtree_metadata.name)
for p in installed.rglob('*'):
    rel=p.relative_to(installed)
    if '__pycache__' in rel.parts or p.suffix in ('.pyc','.pyo'):continue
    target=overlay/'arda'/rel
    if p.is_dir():target.mkdir(exist_ok=True)
    elif p.is_file():
        if '.so' in p.name:target.symlink_to(p)
        else:shutil.copy2(p,target)
for p in baseline.glob('*.libs'):
    if p.is_dir():(overlay/p.name).symlink_to(p)
for name,pin in m['modules'].items():
    target=overlay/'arda'/name
    assert sha(target)==pin['baseline_sha256']
    shutil.copyfile(stage/'source/src/arda'/name,target)
    assert sha(target)==pin['patched_sha256']
assert all(sha(installed/name)==h for name,h in before.items())
files={str(p.relative_to(overlay)):dict(sha256=sha(p),native_symlink=p.is_symlink())
       for p in overlay.rglob('*') if p.is_file()}
vendors={p.name:{str(f.relative_to(p.resolve())):sha(f) for f in p.resolve().rglob('*') if f.is_file()}
         for p in overlay.iterdir() if p.is_symlink() and p.is_dir()}
(stage/'overlay_manifest.json').write_text(json.dumps(dict(source=m,files=files,
    vendor_libraries=vendors,installed_package_unchanged=True),indent=2)+'\n')
PY
