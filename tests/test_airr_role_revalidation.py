"""Retained raw PSM fields and identifier bijections cannot leak evaluation truth."""
import importlib.util
from pathlib import Path
import polars as pl
import pytest


def driver():
    path=Path(__file__).parents[1]/'scripts/revalidate_airr_background_membership.py'
    spec=importlib.util.spec_from_file_location('role_revalidation',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def test_raw_call_projection_strips_truth_and_preserves_acquisition_ids():
    from vdjaqp.readout import _SCHEMA
    m=driver()
    row=dict(peptide_id='solution_trypsin:psm:0',spectrum_id='solution_trypsin:scan=1',
        native_spectrum_reference='scan=1',seq='ACDEFG',local_confidence=[float('nan')]*6,
        search_type='db',alc=float('nan'),alc_rank=1.,n_reads=1,charge=2,rt=1.,mz=500.,hyperscore=10.)
    raw=pl.DataFrame([row],schema=_SCHEMA)
    annotated=raw.with_columns(correct=pl.lit(True),matching_protein_count=pl.lit(12),
                              status=pl.lit('evaluated'),protein_ambiguous=pl.lit(True))
    assert m.raw_calls(annotated,'solution_trypsin').equals(raw)
    assert set(m.raw_calls(annotated,'solution_trypsin').columns)==set(_SCHEMA)
    with pytest.raises(ValueError,match='run_id'):m.raw_calls(annotated,'gel_trypsin')
    with pytest.raises(ValueError,match='missing raw'):m.raw_calls(annotated.drop('seq'),'solution_trypsin')


def test_mapping_preserves_alias_rows_and_role_independent_of_physical_flag():
    m=driver();raw=pl.DataFrame(dict(clone_id=['source_ref','genuine_a','genuine_b'],
        seq=['ACDEFG']*3,nt_seq=['GCT']*3,annotation=['x','y','z']))
    mapping=pl.DataFrame(dict(clone_id=['genuine_b','source_ref','genuine_a'],
        entry_id=['neutral3','neutral1','neutral2'],physical=[True,True,True]))
    result=m.mapped_candidates(raw,mapping,['source_ref'])
    assert result['background_member'].to_list()==[False,True,True]
    assert result['clone_id'].to_list()==['neutral1','neutral2','neutral3']
    assert result.select('seq','nt_seq','annotation').equals(raw.select('seq','nt_seq','annotation'))
    assert result.height==3  # Equal sequence aliases remain separate rows.
    alternate=mapping.with_columns(physical=pl.lit(False))
    assert m.mapped_candidates(raw,alternate,['source_ref']).equals(result)
    with pytest.raises(ValueError,match='bijection'):
        m.mapped_candidates(raw,mapping.with_columns(entry_id=pl.lit('same')),['source_ref'])


def test_all_physical_metrics_do_not_invent_absent_candidate_negatives():
    m=driver();frame=pl.DataFrame(dict(clone_id=['ref','a','b'],score=[3.,2.,1.],confirmed_nats=[25.,0.,0.]))
    rows=m.scoped_metrics(frame,{'ref','a','b'},{'ref'})
    assert rows[0]['metric_scope']=='physical_record_presence'
    assert rows[0]['n_benchmark_negatives']==0 and rows[0]['average_precision'] is None and rows[0]['roc_auc'] is None
    assert rows[1]['metric_scope']=='reference_retrieval' and rows[1]['n_benchmark_negatives']==2
    assert 'physically present' in rows[1]['negative_definition']
    panel=m.scoped_metrics(frame,{'ref'},{'ref'})
    assert len(panel)==1 and panel[0]['average_precision']==1.


def test_original_study_calls_preserve_native_protease_run_ids():
    from vdjaqp.readout import _SCHEMA
    m=driver();rows=[]
    for enzyme in ['Trypsin','Chymotrypsin']:
        rows.append(dict(peptide_id=enzyme+':psm:0',spectrum_id=enzyme+':scan=1',
            native_spectrum_reference='scan=1',seq='ACDEFG',local_confidence=[float('nan')]*6,
            search_type='db',alc=float('nan'),alc_rank=1.,n_reads=1,charge=2,rt=1.,mz=500.,hyperscore=10.))
    frame=pl.DataFrame(rows,schema=_SCHEMA).with_columns(pl.Series('enzyme',['Trypsin','Chymotrypsin']),correct=pl.lit(True))
    actual=m.study_calls(frame)
    assert actual.equals(frame.drop('correct'))
    with pytest.raises(ValueError,match='both named'):m.study_calls(frame.head(1))
