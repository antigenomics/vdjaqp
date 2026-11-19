"""Derive publication macros from final tables/JSON, never from manuscript literals.

Run with a source-root and one or more keyed TSV selector manifests. Table
selection and numerical aggregation use Polars; the loop visits specifications,
not scientific observations. Outputs are requested publication deliverables.
"""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path

import polars as pl


def export(root: Path, manifests: list[Path], output: Path) -> None:
    specs = pl.concat([pl.read_csv(p, separator='\t', infer_schema=False) for p in manifests])
    if specs['key'].null_count() or specs['key'].n_unique() != specs.height:
        raise ValueError('macro keys must be present and unique')
    rows = []
    for spec in specs.iter_rows(named=True):
        source = root / spec['source']
        selector = json.loads(spec['selector'] or '{}')
        aggregation = spec['aggregation']
        if source.suffix == '.json':
            if selector or aggregation not in ('identity', 'count'):
                raise ValueError('JSON sources require an empty selector and identity or count')
            value = json.loads(source.read_text())
            for part in spec['field'].split('.'):
                value = value[int(part)] if isinstance(value, list) else value[part]
            if aggregation == 'count':  # number of entries of a JSON list
                value = len(value)
        else:
            table = pl.read_csv(source, separator='\t', infer_schema=False)
            for field, selected in selector.items():
                table = table.filter(pl.col(field) == str(selected))
            if aggregation == 'count':
                value = table.height
            elif aggregation == 'n_unique':
                value = table[spec['field']].n_unique()
            else:
                values = table[spec['field']].cast(pl.Float64, strict=True)
                if values.null_count() or not values.is_finite().all():
                    raise ValueError(f"missing/nonfinite source: {spec['key']}")
                if aggregation == 'identity':
                    if len(values) != 1:
                        raise ValueError(f"nonunique source: {spec['key']}")
                    value = values.item()
                elif re.fullmatch(r'sorted_index_\d+', aggregation):
                    index = int(aggregation.removeprefix('sorted_index_'))
                    value = values.sort()[index]
                elif aggregation in ('median', 'max', 'min', 'sum', 'mean'):
                    value = getattr(values, aggregation)()
                else:
                    raise ValueError(f'unsupported aggregation: {aggregation}')
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"nonnumeric/nonfinite source: {spec['key']}")
        rows.append((spec['key'], format(value, spec['format']), spec['source'],
                     spec['selector'], spec['field'], aggregation, spec['format'],
                     spec['units'] or '', spec['denominator_note'] or '',
                     hashlib.sha256(source.read_bytes()).hexdigest()))
    output.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(rows, schema=['key', 'value', 'source', 'selector', 'field',
                              'aggregation', 'format', 'units', 'denominator_note',
                              'source_sha256'], orient='row').sort('key').write_csv(output, separator='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, action='append', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    export(args.source_root, args.manifest, args.output)
