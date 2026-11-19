"""Scientific equality must retain nested values/IDs while ignoring row order."""
import importlib.util
from pathlib import Path

import polars as pl


def test_full_row_canonical_comparison():
    path=Path(__file__).parents[1]/'scripts/compare_population_replays.py'
    spec=importlib.util.spec_from_file_location('population_replay_comparison',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    frame=pl.DataFrame(dict(clone_id=['a','b'],seq=['CQQF','CASW'],labels=[[4,6,6,7],[4,5,6,7]],score=[1.5,2.]))
    assert module.canonical(frame).equals(module.canonical(frame.reverse()))
    changed=frame.with_columns(pl.when(pl.col('clone_id')=='a').then(1.6).otherwise(pl.col('score')).alias('score'))
    assert not module.canonical(frame).equals(module.canonical(changed))
    assert not module.canonical(frame).equals(module.canonical(frame.with_columns(pl.col('clone_id')+'x')))


def test_integer_delta_keeps_unit_above_float_precision():
    path=Path(__file__).parents[1]/'scripts/compare_population_replays.py'
    spec=importlib.util.spec_from_file_location('population_replay_counts',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    frame=pl.DataFrame(dict(units=['count','fraction'],fixed1=['9007199254740992','0.5'],
                           fixed2=['9007199254740993','0.5'],alternate=['9007199254740994','0.75']))
    result=module.endpoint_deltas(frame)
    assert result.filter(pl.col('units')=='count').select('same_seed_delta','alternate_minus_fixed1','absolute_difference').row(0)==('1','2','2')
    assert result.filter(pl.col('units')=='fraction')['relative_change'].item()==0.5


def test_schema_bridge_rejects_other_producer_changes(tmp_path):
    import json
    repo=Path(__file__).parents[1]
    path=repo/'scripts/compare_population_replays.py'
    spec=importlib.util.spec_from_file_location('population_replay_source',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    # Invented source pins and exact table evidence keep this test data-free.
    proof=tmp_path/'equivalence.tsv'
    pl.DataFrame([dict(run=run,geometry=mode,table=table,rows=1,schema_and_values_equal=True)
                  for run in ('fixed1a','fixed1b') for mode in ('amplicon','rnaseq')
                  for table in ('catalog','candidates','search_entries')]).write_csv(proof,separator='\t')
    pins=dict(repair='polars-nullable-inference-only', original_sha256='a'*64,
              repaired_sha256='b'*64, evidence_file=proof.name,
              evidence_sha256=module.digest(proof))
    receipt=tmp_path/'source_equivalence.json';receipt.write_text(json.dumps(pins))
    key='scripts/smoke_repertoire.py'
    original={key:pins['original_sha256'],'vdjaqp/bench.py':'same'}
    repaired={**original,key:pins['repaired_sha256']}
    manifests=[{'software_sha256':x} for x in (original,original,repaired)]
    assert module.schema_source_equivalence(receipt,manifests)
    manifests[2]={'software_sha256':{**repaired,'vdjaqp/bench.py':'unexpected'}}
    assert not module.schema_source_equivalence(receipt,manifests)
