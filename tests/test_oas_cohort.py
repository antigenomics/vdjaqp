"""Verify metadata selection and sequence-integrity attrition on synthetic records."""
import json
from pathlib import Path
import subprocess
import sys

import polars as pl


def test_cohort_selection_and_fresh_output(tmp_path):
    metadata = tmp_path / 'metadata.tsv'
    pl.DataFrame({'study': ['Tipton_2015'] * 4, 'vaccine': ['Flu'] * 4,
                  'subject': [f'Subject-FLU{i}' for i in range(1, 5)]}).write_csv(metadata, separator='\t')
    base = dict(sequence='ACDE', nt_sequence='', v_gene='IGHV1', j_gene='IGHJ1',
                length=4, region_ids=[0, 1, 2, 3], study='Tipton_2015',
                donor='Tipton_2015|Subject-FLU1', disease='None')
    changes = [{}, {'disease': 'SARS-COV-2', 'donor': 'other'},
               {'donor': 'Tipton_2015|Subject-SLE1'}, {'sequence': 'ACDX'},
               {'v_gene': 'IGKV1'}, {'length': 3}, {'region_ids': [0, 1, None, 3]}]
    raw = tmp_path / 'raw.parquet'
    pl.DataFrame([dict(base, **change, clone_id=str(i)) for i, change in enumerate(changes)]).write_parquet(raw)
    out = tmp_path / 'selected'
    command = [sys.executable, str(Path(__file__).parents[1] / 'scripts/select_oas_cohort.py'),
               '--input', str(raw), '--influenza-metadata', str(metadata), '--output', str(out)]
    subprocess.run(command, check=True, capture_output=True, text=True)
    selected = pl.read_parquet(out / 'records.parquet')
    assert set(selected['clone_id']) == {'0', '1'}
    assert set(selected['cohort']) == {'covid_disease', 'influenza_vaccination_associated'}
    assert json.loads((out / 'provenance.json').read_text())['records'] == 2
    attrition = pl.read_csv(out / 'attrition.tsv', separator='\t')
    assert attrition.filter(pl.col('stage') == 'metadata_eligible')['records'].sum() == 6
    assert attrition.filter(pl.col('stage') == 'region_labels_present')['records'].sum() == 2
    repeated = subprocess.run(command, capture_output=True, text=True)
    assert repeated.returncode != 0 and 'Output must be new' in repeated.stderr
