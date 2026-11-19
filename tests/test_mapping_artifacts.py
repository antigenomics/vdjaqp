"""Placement artifacts are separate from evidence and preserve native coordinates."""
import hashlib,json
from types import SimpleNamespace
import polars as pl
import pytest
from scripts import workflow_real


def test_native_mapping_artifact_has_offsets_and_separate_substitution_rows(tmp_path):
    clones=pl.DataFrame({'clone_id':['query'],'seq':['AACARDRW']})
    peptides=pl.DataFrame({'peptide_id':['observation'],'seq':['CARERW'],'source_file':['private_source']})
    matches=pl.DataFrame({'clone_id':['query'],'peptide_id':['observation'],'start':[2],'end':[8],'n_subs':[1],'mismatch_positions':[[3]]})
    catalog=tmp_path/'queries.parquet';clones.write_parquet(catalog)
    ev=SimpleNamespace(clonotypes=clones,peptides=peptides,matches=matches)
    out=tmp_path/'mapping'
    workflow_real.export_mappings(ev,out,catalog)
    got=pl.read_parquet(out/'placements.parquet')
    assert got.select('start','end','n_subs').row(0)==(2,8,1)
    assert set(got.columns)=={'placement_id','clone_id','peptide_id','start','end','n_subs'}
    assert pl.read_parquet(out/'substitutions.parquet').row(0)==(0,3)
    assert 'source_file' not in pl.read_parquet(out/'observations.parquet').columns
    receipt=json.loads((out/'mapping_receipt.json').read_text())
    assert receipt['query_catalog_sha256']==hashlib.sha256(catalog.read_bytes()).hexdigest()
    assert receipt['contains_evidence_or_scores'] is False
    assert ev.matches.equals(matches)
    with pytest.raises(FileExistsError):workflow_real.export_mappings(ev,out,catalog)


@pytest.mark.parametrize('n', [0, 1])
def test_empty_or_exact_placements_have_typed_empty_substitutions(tmp_path, n):
    clones=pl.DataFrame({'clone_id':['query'],'seq':['CARDRW']})
    peptides=pl.DataFrame({'peptide_id':['p'],'seq':['CARDRW']})
    matches=pl.DataFrame({'clone_id':['query']*n,'peptide_id':['p']*n,
        'start':[0]*n,'end':[6]*n,'n_subs':[0]*n,'mismatch_positions':[[]]*n},
        schema={'clone_id':pl.String,'peptide_id':pl.String,'start':pl.Int32,
                'end':pl.Int32,'n_subs':pl.Int16,'mismatch_positions':pl.List(pl.Int16)})
    catalog=tmp_path/'catalog.parquet';clones.write_parquet(catalog)
    out=tmp_path/'mapping'
    workflow_real.export_mappings(SimpleNamespace(clonotypes=clones,peptides=peptides,matches=matches),out,catalog)
    assert pl.read_parquet(out/'placements.parquet').height == n
    substitutions=pl.read_parquet(out/'substitutions.parquet')
    assert substitutions.height == 0
    assert substitutions.schema == {'placement_id':pl.UInt32,'peptide_offset':pl.Int16}


@pytest.mark.parametrize('n', [0, 1])
def test_partial_intervals_have_placements_without_evidence(tmp_path, n):
    queries=pl.DataFrame({'clone_id':['partial']*n,'seq':['AACARDRW']*n,
        'identity_eligible':[False]*n},schema={'clone_id':pl.String,'seq':pl.String,'identity_eligible':pl.Boolean})
    catalog=tmp_path/'partial.parquet';queries.write_parquet(catalog)
    peptides=pl.DataFrame({'peptide_id':['p'],'seq':['CARDRW']})
    receipt=workflow_real.export_partial_mappings(catalog,peptides,tmp_path/'placements',1)
    assert receipt['query_role']=='partial_compatibility_only'
    assert receipt['contains_evidence_or_scores'] is False
    assert receipt['placements']==n
    if n:
        assert pl.read_parquet(tmp_path/'placements/placements.parquet').select('start','end').row(0)==(2,8)
    queries.with_columns(pl.lit(True).alias('identity_eligible')).write_parquet(catalog)
    if n:
        with pytest.raises(ValueError,match='identity_eligible'):
            workflow_real.export_partial_mappings(catalog,peptides,tmp_path/'rejected',1)


@pytest.mark.parametrize('eligibility', [[True], [None], [0], ['false']])
def test_partial_query_eligibility_rejected_before_matching(tmp_path, monkeypatch, eligibility):
    from vdjaqp.match import PeptideMatcher
    def forbidden(*args, **kwargs):
        raise AssertionError('invalid catalog reached native matcher')
    monkeypatch.setattr(PeptideMatcher,'from_clonotypes',forbidden)
    queries=pl.DataFrame({'clone_id':['q'],'seq':['CARDRW'],'identity_eligible':eligibility})
    catalog=tmp_path/'partial.parquet';queries.write_parquet(catalog)
    with pytest.raises(ValueError,match='identity_eligible'):
        workflow_real.export_partial_mappings(catalog,pl.DataFrame(),tmp_path/'invalid',1)
    assert not (tmp_path/'invalid').exists()
