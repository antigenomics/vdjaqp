"""Publication selection must fail ambiguous rows and preserve derived values."""
import importlib.util
import json
from pathlib import Path

import polars as pl
import pytest

spec = importlib.util.spec_from_file_location('numbers_export', Path(__file__).parents[1] / 'scripts/export_publication_numbers.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_selected_aggregation_and_rejection(tmp_path):
    pl.DataFrame({'group': ['a', 'a', 'b'], 'value': [1, 3, 9]}).write_csv(tmp_path/'data.tsv', separator='\t')
    manifest = tmp_path/'sources.tsv'
    row = dict(key='median', source='data.tsv', selector='{"group":"a"}', field='value', aggregation='median', format='.1f', units='count', denominator_note='two observations')
    pl.DataFrame([row]).write_csv(manifest, separator='\t')
    module.export(tmp_path, [manifest], tmp_path/'out.tsv')
    assert pl.read_csv(tmp_path/'out.tsv', separator='\t')['value'].item() == 2.0
    row['aggregation'] = 'sorted_index_1'
    pl.DataFrame([row]).write_csv(manifest, separator='\t')
    module.export(tmp_path, [manifest], tmp_path/'index.tsv')
    assert pl.read_csv(tmp_path/'index.tsv', separator='\t')['value'].item() == 3.0
    row['aggregation'] = 'identity'
    pl.DataFrame([row]).write_csv(manifest, separator='\t')
    with pytest.raises(ValueError, match='nonunique'):
        module.export(tmp_path, [manifest], tmp_path/'bad.tsv')

    row.update(selector='{}', field='group', aggregation='n_unique')
    pl.DataFrame([row]).write_csv(manifest, separator='\t')
    module.export(tmp_path, [manifest], tmp_path/'unique.tsv')
    assert pl.read_csv(tmp_path/'unique.tsv', separator='\t')['value'].item() == 2.0
    row['aggregation'] = 'identity'
    (tmp_path/'data.json').write_text(json.dumps({'yield': [17]}))
    row.update(source='data.json', selector='{}', field='yield.0')
    pl.DataFrame([row]).write_csv(manifest, separator='\t')
    module.export(tmp_path, [manifest], tmp_path/'json.tsv')
    assert pl.read_csv(tmp_path/'json.tsv', separator='\t')['value'].item() == 17.0
    with pytest.raises(ValueError, match='unique'):
        module.export(tmp_path, [manifest, manifest], tmp_path/'duplicate.tsv')
