"""Compare three fresh full workflow deliveries without recomputing science.

The fixed pair retains every Parquet field/ID and compares exact row values,
independent of row order. Compressed-container bytes and runtime/path receipts
are not scientific equality criteria; decoded FASTQ hashes are compared instead.
Alternate-seed deltas cover all existing aggregate tables, including changed
gene populations, with no invented stability threshold or missing-row imputation.
"""
import argparse
import hashlib
import json
from pathlib import Path

import polars as pl


def digest(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def canonical(frame):
    """Hash-sort rows, then compare full values; hash equality alone never passes."""
    return frame.with_columns(frame.hash_rows().alias('__row_order')).sort('__row_order').drop('__row_order')


def endpoint_deltas(frame):
    """Subtract counts as integers; floats are used only for fractional changes."""
    parts=[]
    integer_units=pl.col('units').is_in(['count','residues'])
    for integer,dtype in ((True,pl.Int64),(False,pl.Float64)):
        part=frame.filter(integer_units if integer else ~integer_units)
        part=part.with_columns(pl.col(c).cast(dtype,strict=True).alias('_'+c) for c in ('fixed1','fixed2','alternate'))
        part=part.with_columns((pl.col('_fixed2')-pl.col('_fixed1')).alias('same_seed_delta'),
            (pl.col('_alternate')-pl.col('_fixed1')).alias('alternate_minus_fixed1'))
        part=part.with_columns(pl.col('alternate_minus_fixed1').abs().alias('absolute_difference'),
            pl.when(pl.col('_fixed1')!=0).then(pl.col('alternate_minus_fixed1').cast(pl.Float64)/pl.col('_fixed1').cast(pl.Float64)).otherwise(None).alias('relative_change'))
        parts.append(part.drop('_fixed1','_fixed2','_alternate').with_columns(
            pl.col(c).cast(pl.String) for c in ('same_seed_delta','alternate_minus_fixed1','absolute_difference')))
    return pl.concat(parts,how='vertical_relaxed')


def schema_source_equivalence(path,manifests):
    """Accept only the documented, fully tested nullable-constructor repair."""
    receipt=json.loads(path.read_text())
    proof_path=path.parent/receipt['evidence_file']
    if digest(proof_path)!=receipt['evidence_sha256']:raise ValueError('changed schema-equivalence evidence')
    proof=pl.read_csv(proof_path,separator='\t')
    expected={(run,mode,table) for run in ('fixed1a','fixed1b') for mode in ('amplicon','rnaseq') for table in ('catalog','candidates','search_entries')}
    if (proof.height!=12 or proof['schema_and_values_equal'].null_count() or proof['rows'].null_count()
            or set(proof.select('run','geometry','table').iter_rows())!=expected
            or not proof['schema_and_values_equal'].all() or not (proof['rows']>0).all()):
        raise ValueError('complete exact successful-constructor equivalence required')
    a,b,c=[m['software_sha256'] for m in manifests]
    changed={k for k in set(a)|set(c) if a.get(k)!=c.get(k)}
    filename='scripts/smoke_repertoire.py'
    return (receipt['repair']=='polars-nullable-inference-only' and a==b and changed=={filename}
            and a[filename]==receipt['original_sha256'] and c[filename]==receipt['repaired_sha256'])


def aggregates(root):
    keys = dict(component_counts=['component','locus'])
    for mode in ('amplicon','rnaseq'):
        keys.update({mode+'_identity_metrics':['quantity'],
                     mode+'_component_recovery':['component'],
                     mode+'_locus_coverage':['locus','region','evidence_level'],
                     mode+'_component_coverage':['component','locus','region','evidence_level'],
                     mode+'_gene_coverage':['locus','region','gene_call','evidence_level']})
        for enzyme in ('Trypsin','Chymotrypsin'):
            keys[mode+'_'+enzyme+'_readout_metrics']=['quantity']
    dimensions=['component','locus','region','evidence_level','gene_call','quantity']
    frames=[]
    for name,key in keys.items():
        table=pl.read_parquet(root/(name+'.parquet'))
        if table.select(key).unique().height!=table.height:
            raise ValueError('nonunique scientific aggregate keys: '+name)
        numeric=[c for c,t in table.schema.items() if t.is_numeric()]
        for column in numeric:
            unit=('residues' if column in ('covered_residues','residue_denominator') else
                  'fraction' if table.schema[column].is_float() else 'count')
            frame=table.select(*key,pl.col(column).cast(pl.String).alias('value')).with_columns(
                pl.lit(name).alias('table'),pl.lit(column).alias('field'),pl.lit(unit).alias('units'))
            frames.append(frame.with_columns(pl.lit(None,dtype=pl.String).alias(c) for c in dimensions if c not in key))
    classification=pl.read_csv(root/'classification_statistics.tsv',separator='\t',infer_schema_length=0)
    if classification.select('geometry','quantity').unique().height!=36 or classification.height!=36:
        raise ValueError('complete two-geometry classification endpoints required')
    for mode in ('amplicon','rnaseq'):
        frame=classification.filter(pl.col('geometry')==mode).select('quantity','value','units').with_columns(
            pl.lit(mode+'_classification').alias('table'),pl.lit('value').alias('field'))
        frames.append(frame.with_columns(pl.lit(None,dtype=pl.String).alias(c) for c in dimensions if c!='quantity'))
    return pl.concat([f.select('table',*dimensions,'field','units','value') for f in frames],how='vertical_relaxed')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fixed1',type=Path);parser.add_argument('fixed2',type=Path)
    parser.add_argument('alternate',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--schema-equivalence',type=Path,help='reviewed evidence for the sole nullable-constructor repair')
    args=parser.parse_args()
    if args.output.exists():parser.error('fresh output directory required')
    roots=[p/'population' if (p/'population').is_dir() else p for p in (args.fixed1,args.fixed2,args.alternate)]
    manifests=[json.loads((r/'metric_integration.json').read_text()) for r in roots]
    for root,manifest in zip(roots,manifests):
        retained=json.loads((root/'retention.json').read_text())
        if retained.get('removed_artifact_sha256') is None or 'classification_statistics.tsv' not in retained['retained_sha256']:
            raise ValueError('complete final retention receipt required')
        if any(not (root/name).is_file() or digest(root/name)!=pin for name,pin in retained['retained_sha256'].items()):
            raise ValueError('retained final artifact differs from producing receipt')
        expected={name for name in manifest['artifact_sha256'] if name.endswith('.parquet')}
        actual={str(p.relative_to(root)) for p in root.rglob('*.parquet')}
        if not expected or expected!=actual:raise ValueError('producer-declared final scientific inventory mismatch')
        if any(digest(root/name)!=manifest['artifact_sha256'][name] for name in expected):
            raise ValueError('final scientific table differs from its producing receipt')
    files=[{str(p.relative_to(r)) for p in r.rglob('*.parquet')} for r in roots[:2]]
    if files[0]!=files[1]:raise ValueError('fixed runs have different final scientific table sets')
    verdicts=[]
    for name in sorted(files[0]):
        a,b=[pl.read_parquet(r/name) for r in roots[:2]]
        equal=a.schema==b.schema and canonical(a).equals(canonical(b))
        verdicts.append(dict(table=name,rows_fixed1=a.height,rows_fixed2=b.height,
                             exact_values_equal=equal,columns=a.width,
                             file_sha256_fixed1=digest(roots[0]/name),
                             file_sha256_fixed2=digest(roots[1]/name)))
    code_equal=manifests[0]['software_sha256']==manifests[1]['software_sha256']==manifests[2]['software_sha256']
    schema_equivalent=schema_source_equivalence(args.schema_equivalence,manifests) if args.schema_equivalence else False
    configurations=[m['scenario_inputs'] for m in manifests]
    fixed_configuration_equal=configurations[0]==configurations[1]
    nonseed=[{k:v for k,v in c.items() if k not in ('seed','read_seed')} for c in configurations]
    alternate_only_seeds_change=nonseed[0]==nonseed[2] and all(configurations[0][k]!=configurations[2][k] for k in ('seed','read_seed'))
    reads=[json.loads((r/'read_run_manifest.json').read_text()) for r in roots]
    logical=[{x['scenario']:x['logical_fastq'] for x in run} for run in reads]
    native=[]
    for root in roots:
        producer=json.loads((root/'workflow_producer.json').read_text())
        pins={k:v for k,v in producer.items() if k!='read_error_seed'}
        for mode in ('amplicon','rnaseq'):
            for enzyme in ('Trypsin','Chymotrypsin'):
                name=mode+'_'+enzyme
                summary=json.loads((root/(name+'_native_summary.json')).read_text())
                pins[name]={k:summary[k] for k in ('selected_mssimulator_sha256','search_export_executables','settings','threads')}
        native.append(pins)
    native_equal=native[0]==native[1]==native[2]
    classification=[pl.read_csv(r/'classification_statistics.tsv',separator='\t',infer_schema_length=0) for r in roots[:2]]
    classification_equal=classification[0].schema==classification[1].schema and canonical(classification[0]).equals(canonical(classification[1]))
    key=['table','component','locus','region','evidence_level','gene_call','quantity','field','units']
    merged=None
    for label,root in zip(('fixed1','fixed2','alternate'),roots):
        frame=aggregates(root).rename({'value':label}).with_columns(pl.lit(True).alias(label+'_row_present'))
        merged=frame if merged is None else merged.join(frame,on=key,how='full',coalesce=True,nulls_equal=True,validate='1:1')
    merged=endpoint_deltas(merged).sort(key)
    args.output.mkdir(parents=True)
    pl.DataFrame(verdicts).write_csv(args.output/'same_seed_tables.tsv',separator='\t')
    merged.write_csv(args.output/'endpoint_deltas.tsv',separator='\t')
    pl.concat([aggregates(r).with_columns(pl.lit(label).alias('run')) for label,r in zip(('fixed1','fixed2','alternate'),roots)]).write_csv(args.output/'scientific_statistics.tsv',separator='\t')
    report=dict(scientific_tables=len(verdicts),same_seed_exact_all_tables=all(r['exact_values_equal'] for r in verdicts),
        same_seed_classification_tsv_exact=classification_equal,
        same_seed_logical_fastq_equal=logical[0]==logical[1],alternate_logical_fastq_changed=logical[0]!=logical[2],
        producing_software_hashes_equal=code_equal,endpoint_rows=merged.height,
        reviewed_schema_only_source_equivalence=schema_equivalent,
        schema_equivalence_receipt_sha256=digest(args.schema_equivalence) if args.schema_equivalence else None,
        native_runtime_and_settings_equal=native_equal,native_runtime_pins=native,
        fixed_configuration_equal=fixed_configuration_equal,alternate_only_seeds_change=alternate_only_seeds_change,
        scenario_inputs=configurations,
        read_configurations=[[dict(geometry=x['scenario'],seed=x['seed'],cpus=x['cpus'],pairs=x['pairs']) for x in run] for run in reads],
        logical_fastq=logical,missing_endpoint_policy='outer join; absent populations stay missing; no zero imputation',
        relative_change='(alternate-fixed1)/fixed1; undefined if fixed1 zero or either observation unavailable',
        interpretation='one alternate seed is a bounded stability observation, not an uncertainty estimate; no tolerance selected')
    (args.output/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if isinstance(v,(int,bool))}))
    if not (report['same_seed_exact_all_tables'] and classification_equal and report['same_seed_logical_fastq_equal'] and (code_equal or schema_equivalent) and native_equal and fixed_configuration_equal and alternate_only_seeds_change):raise SystemExit(1)


if __name__=='__main__':main()
