"""Qualify the native OAS aggregation on known synthetic partitions, not scientific data.

Usage: python tests/check_oas_digest.py BINARY NEW_OUTPUT_DIRECTORY [CPU_COUNT]
"""
import collections
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time


def check(binary, output, cpus=4):
    output.mkdir(parents=True, exist_ok=False)
    proteins = ('AKPRAK', 'AFPLWYA', 'DADD')
    partitions = {
        'Trypsin': (('AKPR', 'AK'), ('AFPLWYA',), ('DADD',)),
        'Chymotrypsin': (('AKPRAK',), ('AFPL', 'W', 'Y', 'A'), ('DADD',)),
        'Asp-N': (('AKPRAK',), ('AFPLWYA',), ('DA', 'D', 'D')),
    }
    repeats = 10000
    source = output / 'synthetic.tsv'
    with source.open('w') as f:
        for _ in range(repeats):
            for i, protein in enumerate(proteins):
                # Middle protein deliberately has no inferred regions; native digestion must retain it.
                bounds = [0, len(protein) if i != 1 else 0] + [0] * 12
                f.write('\t'.join(map(str, [i % 2, i, protein, *bounds])) + '\n')
    elapsed = {}
    for threads in sorted({1, cpus}):
        start = time.monotonic()
        subprocess.run([str(binary), str(source), str(output / f't{threads}-'), str(threads)], check=True)
        elapsed[str(threads)] = time.monotonic() - start
    for suffix in ('records', 'lengths', 'native', 'aligned'):
        assert (output / f't1-{suffix}.tsv').read_bytes() == (output / f't{cpus}-{suffix}.tsv').read_bytes()
    def rows(suffix):
        with (output / f't1-{suffix}.tsv').open() as f:
            return list(csv.DictReader(f, delimiter='\t'))
    assert {int(r['cohort_id']): int(r['protein_records']) for r in rows('records')} == {0: 2*repeats, 1: repeats}
    expected = collections.Counter()
    native, aligned = {}, {}
    for i, protein in enumerate(proteins):
        for c in range(1, len(protein)):
            key = (i % 2, len(protein)//10, (c-1)//5)
            native.setdefault(key, collections.Counter())['available_bonds'] += repeats
            if i != 1:
                aligned.setdefault((i % 2, i, 'fwr1', c), collections.Counter())['available_bonds'] += repeats
        for enzyme, parts in partitions.items():
            fragments = parts[i]
            for width in range(1, min(4, len(fragments))+1):
                for start in range(len(fragments)-width+1):
                    length = sum(map(len, fragments[start:start+width]))
                    expected[(i % 2, enzyme, width-1, length)] += repeats
            c = 0
            for fragment in fragments[:-1]:
                c += len(fragment)
                native[(i % 2, len(protein)//10, (c-1)//5)][enzyme] += repeats
                if i != 1:
                    aligned[(i % 2, i, 'fwr1', c)][enzyme] += repeats
    actual = {(int(r['cohort_id']), r['enzyme'], int(r['missed_cleavages']), int(r['length_aa'])): int(r['peptide_occurrences']) for r in rows('lengths')}
    assert actual == expected, (actual, expected)
    for suffix, expected_map in [('native', native), ('aligned', aligned)]:
        actual_map = {}
        for row in rows(suffix):
            key = ((int(row['cohort_id']), int(row['length_bin10']), int(row['bond_bin5'])) if suffix == 'native'
                   else (int(row['cohort_id']), int(row['gene_id']), row['region'], int(row['bond_offset'])))
            actual_map[key] = collections.Counter({name: int(row[name]) for name in ['available_bonds', *partitions]})
        assert actual_map == expected_map, suffix
    # Reject malformed geometry rather than treating it as absent annotation.
    bad = output / 'invalid.tsv'
    bad.write_text('0\t0\tAK\t1\t0' + '\t0'*12 + '\n')
    result = subprocess.run([str(binary), str(bad), str(output/'invalid-'), '1'], capture_output=True, text=True)
    assert result.returncode != 0 and 'Incomplete inferred region' in result.stderr
    report = dict(status='complete', records=3*repeats, exact_serial_parallel_equality=True,
                  known_partition_counts_verified=True, unresolved_regions_retained=True,
                  seconds=elapsed, binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                  fixture_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (output/'receipt.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    check(Path(sys.argv[1]).resolve(), Path(sys.argv[2]), int(sys.argv[3]) if len(sys.argv)>3 else 4)
