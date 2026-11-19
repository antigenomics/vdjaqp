"""Rebuild both tiny published-tool scenarios from raw references, without resume.

Purpose: software integration verification, not scientific performance evaluation.
Requires the separately installed, pinned runtimes in envs/.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
from urllib.request import urlopen
import xml.etree.ElementTree as ET


def prepare_references(repo, root, *, all_loci=False, source_dir=None):
    """Fetch only checksum-pinned raw records for a fresh fixture."""
    refs = root / 'references'
    refs.mkdir(parents=True)
    names = ['references.json'] + (['lambda_references.json'] if all_loci else [])
    for filename in names:
        manifest = json.loads((repo / 'envs/generation' / filename).read_text())
        for name, entry in manifest['references'].items():
            response = (Path(source_dir)/name).open('rb') if source_dir is not None else urlopen(entry['url'], timeout=30)
            with response:
                content = response.read(2_000_001)
            if len(content) > 2_000_000 or hashlib.sha256(content).hexdigest() != entry['sha256']:
                raise ValueError(f'raw reference size/hash mismatch: {name}')
            (refs / name).write_bytes(content)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='new directory; must not exist')
    parser.add_argument('--cpus', type=int, default=None)
    parser.add_argument('--reference-dir',type=Path,help='pre-staged raw records; every checksum is verified')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output already exists; choose a new directory')
    repo = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo))
    import vdjaqp
    cpus = vdjaqp.resolve_cpu_budget(args.cpus)
    root = args.output.resolve()
    prepare_references(repo, root, source_dir=args.reference_dir)
    generation = repo / '.tool-envs/generation/bin/python'
    stages = []

    def run(label, interpreter, script, *arguments):
        start = time.perf_counter()
        command = [str(interpreter), str(repo / 'scripts' / script), *map(str, arguments)]
        with (root / (label + '.log')).open('w') as log:
            subprocess.run(command, cwd=repo, stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=300)
        stages.append({'stage': label, 'seconds': time.perf_counter()-start, 'command': command})

    # Generation permits raw references only; stage logs are created after it starts.
    generation_command = [str(generation), str(repo/'scripts/smoke_generate.py'),
                          '--output', str(root)]
    start = time.perf_counter()
    subprocess.run(generation_command, cwd=repo, check=True, timeout=300)
    stages.append({'stage': 'generation', 'seconds': time.perf_counter()-start, 'command': generation_command})
    run('reads', generation, 'smoke_reads.py', '--root', root, '--cpus', cpus)
    summaries = {}
    for scenario in ('amplicon', 'rnaseq'):
        run(scenario+'-repertoire', sys.executable, 'smoke_repertoire.py', scenario, '--root', root, '--threads', cpus)
        from smoke_metrics import load_observed
        bundle=load_observed(root/(scenario+'_observed'),scenario)
        database=root/(scenario+'_observed.fasta')
        database.write_text(''.join('>'+r['entry_id']+'\n'+r['seq']+'\n' for r in bundle['search_entries'].iter_rows(named=True)))
        native = root / ('proteomics-'+scenario)
        run(scenario+'-proteomics', sys.executable, 'smoke_proteomics.py', root/'proteins.fasta',
            root/(scenario+'_observed.fasta'), native, '--threads', cpus, '--keep-spectra')
        scores = root / (scenario+'_scores.parquet')
        run(scenario+'-score', sys.executable, 'smoke_score.py', root/(scenario+'_observed/candidates.parquet'),
            native/'psms.idXML', scores, '--run-id', scenario, '--threads', cpus)
        import polars as pl
        peptides = pl.read_parquet(root/(scenario+'_scores_peptides.parquet'))
        spectra = native/'spectra.mzML'
        identities = {node.attrib['id'] for _, node in ET.iterparse(spectra, events=('end',))
                      if node.tag.endswith('}spectrum')}
        assert set(peptides['native_spectrum_reference']) <= identities
        truth_rows = json.loads((native/'ms2_truth.json').read_text())
        truth_ids = {row['native_spectrum_reference'] for row in truth_rows}
        assert set(peptides['native_spectrum_reference']) <= truth_ids
        assert len(truth_ids) == len(truth_rows)
        report = json.loads(scores.with_suffix('.json').read_text())
        report.update(native_spectrum_ids_verified=True, native_ms2_truth_verified=True,
                      native_ms2_truth_groups=len(truth_rows),
                      artifact_sha256={str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in (root/'truth.json', root/'proteins.fasta',
                              root/(scenario+'_observed.fasta'), native/'psms.idXML',
                              native/'truth.idXML', native/'ms2_truth.json', scores,
                              root/(scenario+'_scores_peptides.parquet'))},
                      discarded_spectra_sha256=hashlib.sha256(spectra.read_bytes()).hexdigest(),
                      discarded_spectra_bytes=spectra.stat().st_size, transient_spectra_removed=True)
        spectra.unlink()
        scores.with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
        summaries[scenario] = report
    (root/'end_to_end.json').write_text(json.dumps({'cpus':cpus, 'stages':stages, 'scenarios':summaries, 'software_sha256': {str(path.relative_to(repo)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted([*repo.glob('vdjaqp/*.py'), *repo.glob('scripts/smoke*.py'), *repo.glob('envs/*/*'), repo/'uv.lock']) if path.is_file()}}, indent=2)+'\n')
    print(json.dumps(summaries, indent=2))


if __name__ == '__main__':
    main()
