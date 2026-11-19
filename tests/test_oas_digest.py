"""Coordinate/provenance checks for native OAS digestion; no scientific recomputation."""
import json
from pathlib import Path
import sys

import polars as pl
import pytest
from scripts import digest_oas_cohort as driver


def inputs():
    proteins = pl.DataFrame({'clone_id': ['a', 'b'], 'cohort': ['test', 'test'], 'sequence': ['AKPR', 'AAAA']})
    rows = []
    for ident, protein in [('a', 'AKPR'), ('b', 'AAAA')]:
        row = dict(sequence_id=ident, sequence=protein, v_call='IGHV1' if ident == 'a' else '')
        for r in driver.REGIONS:
            row.update({r+'_start': '', r+'_end': '', r+'_aa': ''})
        if ident == 'a':
            row.update(fwr1_start='1', fwr1_end='2', fwr1_aa='AK', cdr1_start='3', cdr1_end='4', cdr1_aa='PR')
        rows.append(row)
    return proteins, pl.DataFrame(rows)


def test_coordinates_and_unresolved():
    proteins, annotations = inputs()
    table, cohorts, genes, quality = driver.prepare(proteins, annotations)
    known = table.filter(pl.col('sequence') == 'AKPR').row(0, named=True)
    assert (known['fwr1_start0'], known['fwr1_end0'], known['cdr1_start0'], known['cdr1_end0']) == (0, 2, 2, 4)
    assert table.height == 2 and genes['v_call'].to_list() == ['', 'IGHV1']
    assert quality['defined_regions'].to_list() == [0, 2]
    for change in [dict(fwr1_aa='AA'), dict(cdr1_start='2', cdr1_aa='KPR'), dict(fwr1_end='')]:
        bad = annotations.to_dicts(); bad[0].update(change)
        with pytest.raises(ValueError, match='Inconsistent native'):
            driver.prepare(proteins, pl.DataFrame(bad))
    with pytest.raises(ValueError, match='missing or duplicated'):
        driver.prepare(proteins, annotations.head(1))


def test_native_integer_coordinates_match_airr_text():
    proteins, annotations = inputs()
    coordinates = [r+s for r in driver.REGIONS for s in ('_start', '_end')]
    native = annotations.with_columns(pl.col(coordinates).replace('', None).cast(pl.Int64))
    expected = driver.prepare(proteins, annotations)
    observed = driver.prepare(proteins, native)
    assert all(a.equals(b) for a, b in zip(expected, observed))


def test_driver_retains_all_records_and_integer_totals(tmp_path, monkeypatch):
    proteins, annotations = inputs()
    root = tmp_path/'annotations'; selection = root/'selection'; selection.mkdir(parents=True)
    proteins.write_parquet(selection/'records.parquet'); annotations.write_parquet(root/'annotations.parquet')
    (selection/'provenance.json').write_text(json.dumps({'output_sha256': driver.sha(selection/'records.parquet')}))
    (root/'receipt.json').write_text(json.dumps({'status': 'complete', 'annotations_sha256': driver.sha(root/'annotations.parquet')}))
    binary = tmp_path/'binary'; binary.write_text('explicit test stub')
    prefix = tmp_path/'openms'
    for name in ['lib/libOpenMS.so', 'share/OpenMS/CHEMISTRY/Enzymes.xml']:
        p = prefix/name; p.parent.mkdir(parents=True, exist_ok=True); p.write_text('fixture')
    native = tmp_path/'native.json'
    native.write_text(json.dumps(dict(source_sha256=driver.sha(Path(driver.__file__).with_name('oas_digest.cpp')),
        binary_sha256=driver.sha(binary), library_sha256=driver.sha(prefix/'lib/libOpenMS.so'),
        enzyme_definitions_sha256=driver.sha(prefix/'share/OpenMS/CHEMISTRY/Enzymes.xml'))))
    output = tmp_path/'result'
    monkeypatch.setattr(sys, 'argv', ['digest', '--annotations', str(root), '--binary', str(binary),
        '--native-receipt', str(native), '--openms-prefix', str(prefix), '--output', str(output), '--cpus', '1'])
    monkeypatch.setattr(driver, 'resolve_cpu_budget', lambda _: 1)
    def run(command, **kwargs):
        p = command[2]
        Path(p+'native.tsv').write_text('cohort_id\tlength_bin10\tbond_bin5\tavailable_bonds\tTrypsin\tChymotrypsin\tAsp-N\n0\t0\t0\t6\t0\t0\t0\n')
        Path(p+'aligned.tsv').write_text('cohort_id\tgene_id\tregion\tbond_offset\tavailable_bonds\tTrypsin\tChymotrypsin\tAsp-N\n0\t1\tfwr1\t1\t1\t0\t0\t0\n')
        Path(p+'records.tsv').write_text('cohort_id\tprotein_records\n0\t2\n')
        Path(p+'lengths.tsv').write_text('cohort_id\tenzyme\tmissed_cleavages\tlength_aa\tpeptide_occurrences\n'+''.join(f'0\t{e}\t0\t4\t2\n' for e in driver.ENZYMES))
    monkeypatch.setattr(driver.subprocess, 'run', run)
    driver.main()
    receipt = json.loads((output/'receipt.json').read_text())
    assert receipt['records'] == 2 and receipt['status'] == 'complete'
    assert not (output/'proteins.tsv').exists()
    assert pl.read_csv(output/'length_eligibility.tsv', separator='\t')['length_eligible_occurrences'].sum() == 0
    assert all(driver.sha(output/p) == h for p, h in receipt['output_sha256'].items())
