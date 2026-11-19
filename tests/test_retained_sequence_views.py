"""Sequence displays consume the fresh DAG's artifacts without another model run."""
import hashlib
import json
from types import SimpleNamespace
import polars as pl
import pytest
from scripts.workflow_real import export_mappings
from scripts.export_sequence_evidence import export_retained


def fixture(tmp_path, empty=False):
    inputs=tmp_path/'inputs';(inputs/'ground_truth').mkdir(parents=True)
    targets=inputs/'ground_truth/target_clonotypes.csv'
    pl.DataFrame({'clone_id':['control']}).write_csv(targets)
    clones=pl.DataFrame(dict(clone_id=['control','rna'],seq=['AACILDEWGG','AACALDEWGG'],
        region_ids=[[0,0,4,5,6,6,5,7,7,7]]*2,cdr3_start=[2]*2,cdr3_end=[8]*2,
        sequence_origin=['supplied_control','observed_native_interval']))
    catalog=inputs/'catalog.parquet';clones.write_parquet(catalog)
    context=tmp_path/'context';dest=context/'analysis/gel';dest.mkdir(parents=True)
    clones.write_parquet(dest/'candidates.parquet')
    scores=pl.DataFrame(dict(clone_id=['control','rna'],confirmed_nats=[0.,0.] if empty else [24.,10.],
        junction_nats=[0.,0.] if empty else [12.,4.],framework_only=[False,False]))
    scores.write_parquet(dest/'alignment_union_scores.parquet')
    peptides=pl.DataFrame({'peptide_id':['p'],'seq':['CLLDEW']})
    matches=pl.DataFrame(dict(clone_id=['control','rna'],peptide_id=['p']*2,
        start=[2,2],end=[8,8],n_subs=[0,1],mismatch_positions=[[],[1]]),
        schema_overrides={'mismatch_positions':pl.List(pl.Int16)})
    if empty:matches=matches.head(0)
    mapping=context/'peptide_mapping'
    export_mappings(SimpleNamespace(clonotypes=clones,peptides=peptides,matches=matches),mapping,catalog)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    provenance=dict(status='complete',cases=[dict(case='gel',mapping_receipt_sha256=sha(mapping/'mapping_receipt.json'))],
        raw_sha256={'ground_truth/target_clonotypes.csv':sha(targets)},
        artifact_sha256={str(p.relative_to(context/'analysis')):sha(p) for p in dest.iterdir()})
    (context/'analysis/provenance.json').write_text(json.dumps(provenance))
    display=tmp_path/'display.tsv'
    pl.DataFrame(dict(clone_id=['rna','control'],display_id=['ExampleB','ExampleA'])).write_csv(display,separator='\t')
    return inputs,context,display,tmp_path/'views'


@pytest.mark.parametrize('empty',[False,True])
def test_retained_views_preserve_scores_coordinates_and_full_universe(tmp_path,monkeypatch,empty):
    from vdjaqp.pipeline import Evidence
    def forbidden(*args,**kwargs):raise AssertionError('display attempted scientific recomputation')
    for name in ['build','from_frames','score']:monkeypatch.setattr(Evidence,name,forbidden)
    inputs,context,display,out=fixture(tmp_path,empty)
    export_retained(inputs,context,display,out,2,cpus=1)
    table=pl.read_csv(out/'clonotypes.tsv',separator='\t')
    assert table['display_id'].to_list()==['ExampleB','ExampleA']
    assert table['covered_aa'].to_list()==([0,0] if empty else [5,6])
    assert table['confirmed_nats'].to_list()==([0.,0.] if empty else [10.,24.])
    hits=pl.read_parquet(out/'placements.parquet')
    assert hits.height==(0 if empty else 2)
    if not empty:
        assert hits['start'].to_list()==[2,2]
        assert hits['compatible_candidate_records'].to_list()==[2,2]
    residue=pl.read_csv(out/'residues.tsv',separator='\t').filter(pl.col('display_id')=='ExampleA')
    assert residue.filter(pl.col('in_junction'))['position'].to_list()==list(range(2,8))
    receipt=json.loads((out/'manifest.json').read_text())
    assert receipt['search_or_scoring_repeated'] is False
    assert receipt['candidate_count']==2
    with pytest.raises(FileExistsError):export_retained(inputs,context,display,out,2,cpus=1)


def test_retained_view_rejects_artifact_drift_before_writing(tmp_path):
    inputs,context,display,out=fixture(tmp_path)
    (context/'peptide_mapping/placements.parquet').write_bytes(b'changed')
    with pytest.raises(ValueError,match='artifact drift'):
        export_retained(inputs,context,display,out,2,cpus=1)
    assert not out.exists()


def test_display_selection_change_during_read_is_not_attested_as_used(tmp_path,monkeypatch):
    inputs,context,display,out=fixture(tmp_path)
    original=pl.read_csv
    def changing(path,*args,**kwargs):
        frame=original(path,*args,**kwargs)
        if path==display:
            display.write_text('clone_id\tdisplay_id\ncontrol\tChanged\n')
        return frame
    monkeypatch.setattr(pl,'read_csv',changing)
    with pytest.raises(ValueError,match='input drift'):
        export_retained(inputs,context,display,out,2,cpus=1)
    assert not out.exists()
