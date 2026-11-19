"""Build and qualify the OAS batch wrapper against an existing Linux OpenMS 3.0 prefix.

The installed OpenMS environment is read-only; no packages or protease rules change.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from vdjaqp.resources import resolve_cpu_budget
from scripts.annotate_oas_cohort import sha
from scripts.build_openms_digest_fix import build_flags


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--openms-prefix', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpus', default='all', help='CPU allocation for synthetic serial/parallel qualification')
    args = parser.parse_args()
    cpus = resolve_cpu_budget(args.cpus)
    if args.output.exists():
        parser.error('Build output must be new.')
    prefix = args.openms_prefix.resolve()
    library = prefix/'lib/libOpenMS.so'
    enzymes = prefix/'share/OpenMS/CHEMISTRY/Enzymes.xml'
    if not library.is_file():
        parser.error('This wrapper requires the qualified Linux OpenMS 3.0 environment.')
    source = Path(__file__).with_name('oas_digest.cpp').resolve()
    fixture = source.parents[1]/'tests/check_oas_digest.py'
    pins = {str(p): sha(p) for p in (library, enzymes, source, fixture)}
    _, flags, link, _ = build_flags(prefix)
    args.output.mkdir(parents=True)
    binary = (args.output/'oas_digest').resolve()
    command = flags+['-fopenmp', str(source)]+link+['-o', str(binary)]
    subprocess.run(command, check=True)
    check = args.output/'qualification'
    subprocess.run([sys.executable, str(fixture), str(binary), str(check), str(cpus)], check=True)
    if any(sha(Path(p)) != expected for p, expected in pins.items()):
        raise ValueError('Native source or runtime changed during qualification.')
    receipt = dict(status='complete', command=command, source_sha256=sha(source), binary_sha256=sha(binary),
                   library_sha256=sha(library), enzyme_definitions_sha256=sha(enzymes),
                   qualification=json.loads((check/'receipt.json').read_text()), source_input_sha256=pins)
    (args.output/'native_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    # Synthetic tables can be regenerated; retain the compact exact-check receipt.
    for table in check.glob('*.tsv'):
        table.unlink()
    print(args.output/'native_receipt.json')


if __name__ == '__main__':
    main()
