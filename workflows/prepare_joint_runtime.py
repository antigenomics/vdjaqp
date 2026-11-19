"""Retain compact task-owned qualification proof and raw refs, remove generated bulk.

Run only in the new tracked reproduction runtime, after its baseline jobs end.
Original supplied inputs and pinned environments are outside the deletion roots.
"""
from pathlib import Path
import hashlib
import json
import shutil
import sys

root=Path(sys.argv[1]).resolve()
if root.name!='2026-igprotseq-reproduction-20261112':raise ValueError('unexpected runtime')
proof=root/'qualification';proof.mkdir(exist_ok=False)
old=root/'null/population'
shutil.copytree(old/'references',root/'raw-reference-inputs')
profile=proof/'profile24';profile.mkdir()
for p in old.iterdir():
    if p.is_file() and p.suffix in ('.parquet','.tsv','.json'):shutil.copy2(p,profile/p.name)
for folder in ('amplicon_observed','rnaseq_observed'):
    shutil.copytree(old/folder,profile/folder)
failed=root/'code/tmp/nextflow-work/88/9cccfd5b9fc5c19ba687d947e847df'
shutil.copy2(failed/'population-resources.json',proof/'unpatched_full_failure.json')
truth=json.loads((failed/'population/truth.json').read_text())
chains=truth['chains']
(proof/'unpatched_full_substrate.json').write_text(json.dumps(dict(
    seed=truth['seed'],chain_count=len(chains),antibody_count=len({c['antibody_id'] for c in chains}),
    pairs_per_geometry=sum(c['rna_pairs'] for c in chains),
    truth_sha256=hashlib.sha256((failed/'population/truth.json').read_bytes()).hexdigest(),
    sources=truth.get('reference_sha256'),status='MS stage timed out; no complete scientific recovery output'),indent=2)+'\n')
for name in ('read_run_manifest.json',):shutil.copy2(failed/'population'/name,proof/name)
for mode in ('amplicon','rnaseq'):
    shutil.copy2(failed/'population'/(mode+'_observed')/'receipt.json',proof/(mode+'_native_receipt.json'))
removed=[]
for path in (failed/'population',old):
    removed.append(str(path.relative_to(root)))
    shutil.rmtree(path)
(proof/'cleanup.json').write_text(json.dumps(dict(removed=removed,
    preserved='original-inputs, raw-reference-inputs, pinned environments, final real deliverables and compact qualification proof'),indent=2)+'\n')
