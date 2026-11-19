"""Recount the complete experimental cohort from raw supplied inputs.

The author requested unit count placeholders. Original Parquet inputs remain
unchanged; derived candidates and capture tables label those counts explicitly.
Every invocation computes every context afresh, then replaces only its own
previous completed output. Failed runs never replace the previous result.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

repo = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo))
from vdjaqp import resolve_cpu_budget


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpus', default='all')
    args = parser.parse_args()
    budget = resolve_cpu_budget(args.cpus)
    inputs = repo / 'private_data/igseq_datasets'
    output = repo / 'private_data/publication_analysis'
    marker = '.publication-analysis'
    if output.is_symlink() or (output.exists() and not (output / marker).is_file()):
        parser.error('refusing to replace an unowned or symlink output directory')
    scratch = repo / '.build'
    scratch.mkdir(exist_ok=True)
    with TemporaryDirectory(prefix='analyze-', dir=scratch) as temporary:
        stage = Path(temporary)
        command = [sys.executable, str(repo / 'scripts/profile_build.py'),
                   '--output', str(stage / 'resource.json'), '--project-root', str(repo),
                   '--workload-root', str(stage), '--timeout', '1650',
                   '--rss-limit', '24000000000', '--disk-limit', '10000000000',
                   '--cpus', str(budget), '--', sys.executable,
                   str(repo / 'scripts/workflow_real.py'), '--input-root', str(inputs),
                   '--output', str(stage / 'analysis'), '--cpus', str(budget),
                   '--placeholder-count-one', '--publication']
        result = subprocess.run(command, cwd=repo, timeout=1720)
        if result.returncode:
            for suffix in ('json', 'log'):
                report = stage / f'resource.{suffix}'
                if report.is_file():
                    shutil.copy2(report, repo / f'private_data/publication_analysis_failure.{suffix}')
            raise SystemExit(result.returncode)
        report = json.loads((stage / 'analysis/provenance.json').read_text())
        if report['status'] != 'complete' or len(report['cases']) != 10:
            raise ValueError('experimental recount incomplete; previous outputs preserved')
        for suffix in ('json', 'log'):
            shutil.copy2(stage / f'resource.{suffix}', stage / f'analysis/resource.{suffix}')
        (stage / 'analysis' / marker).write_text('Fresh ten-context experimental publication analysis\n')
        if output.exists():
            output.rename(stage / 'previous')
        try:
            (stage / 'analysis').rename(output)
        except BaseException:
            if (stage / 'previous').exists():
                (stage / 'previous').rename(output)
            raise
    print(f'Fresh experimental publication outputs: {output}', flush=True)


if __name__ == '__main__':
    main()
