"""Exercise the fresh-run driver with an explicit native-CLI stub, never biological inference."""
import json
import sys

import arda
import polars as pl
import pytest

from scripts import annotate_oas_cohort as driver


@pytest.mark.parametrize('omit_record', [False, True])
def test_native_records_preserved_or_failure(tmp_path, monkeypatch, omit_record):
    raw = tmp_path / 'raw.parquet'
    pl.DataFrame([dict(clone_id=str(i), sequence='ACDE', nt_sequence='', v_gene='IGHV1',
                      j_gene='IGHJ1', length=4, region_ids=[0, 1, 2, 3], study='Study',
                      donor='donor', disease='SARS-COV-2') for i in range(2)]).write_parquet(raw)
    metadata = tmp_path / 'metadata.tsv'
    pl.DataFrame({'study': ['Tipton_2015'] * 4, 'vaccine': ['Flu'] * 4,
                  'subject': [f'Subject-FLU{i}' for i in range(1, 5)]}).write_csv(metadata, separator='\t')
    output = tmp_path / 'result'
    monkeypatch.setattr(sys, 'argv', ['annotation', '--input', str(raw), '--influenza-metadata',
                                    str(metadata), '--output', str(output), '--cpus', '1'])
    monkeypatch.setattr(driver, 'resolve_cpu_budget', lambda value: 1)
    monkeypatch.setattr(arda, '__version__', '2.38.7')
    monkeypatch.setattr(arda, 'reference_index_identity', lambda *a, **k: {'fixture': 'reference'})
    run = driver.subprocess.run

    def native(command, **kwargs):
        if 'annotate' not in command:
            return run(command, **kwargs)
        fasta = output / 'proteins.fasta'
        assert fasta.read_text() == '>0\nACDE\n>1\nACDE\n'
        (output / 'annotations.airr.tsv').write_text(
            'sequence_id\tsequence\tlocus\n0\tACDE\t\n' + ('' if omit_record else '1\tACDE\tIGH\n'))

    monkeypatch.setattr(driver.subprocess, 'run', native)
    if omit_record:
        with pytest.raises(ValueError, match='lost or duplicated'):
            driver.main()
        assert not (output / 'receipt.json').exists()
    else:
        driver.main()
        records = pl.read_parquet(output / 'annotations.parquet')
        assert records['locus'].to_list() == ['', 'IGH']  # unresolved remains present
        assert json.loads((output / 'receipt.json').read_text())['records'] == 2
        assert not (output / 'proteins.fasta').exists()
        assert not (output / 'annotations.airr.tsv').exists()
