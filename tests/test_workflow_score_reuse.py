"""Each fresh context computes common scorer state and the full layer once."""
import sys
import json
import numpy as np
import polars as pl
from scripts import workflow_real
from vdjaqp.pipeline import Evidence
from vdjaqp.errormodel import ErrorModel,mass_prior
from vdjaqp import background


def test_publication_reuses_current_run_full_score_without_changing_results(tmp_path,monkeypatch):
    root=tmp_path/'inputs';(root/'ground_truth').mkdir(parents=True)
    ids=[f'Ab_{i:02d}' for i in range(12)]+['rna']
    pl.DataFrame(dict(clone_id=ids[:12],has_single_run=['yes']*8+['no']*4,
        ms_stem=[f'run{i}' for i in range(12)])).write_csv(root/'ground_truth/target_clonotypes.csv')
    clones=pl.DataFrame(dict(clone_id=ids,seq=['CARGDW']*13,region_ids=[[6]*6]*13,
        cdr3_seq=['CARGDW']*13,cdr3_start=[0]*13,cdr3_end=[6]*13,
        v_name=['IGHV1-2']*13,j_name=['IGHJ4']*13,c_name=['IGHG1']*12+['IGHM'],
        background_member=[False]*12+[True]))
    peptides=pl.DataFrame(dict(peptide_id=['run/1000/2/500/10/1/CARGDW'],seq=['CARGDW'],
        local_confidence=[[.99]*6],search_type=['db'],alc=[.99],alc_rank=[1.],n_reads=[1]))
    monkeypatch.setattr(background,'canonical_v',lambda:{})
    monkeypatch.setattr(background,'canonical_j',lambda:{})
    model=ErrorModel(np.array([0.,1.]),np.array([.1,.01]),mass_prior(),0)
    ev=Evidence.from_frames(clones,peptides,error_model=model)
    calls=[];original=ev.background.information
    def information():calls.append(True);return original()
    ev.background.information=information
    monkeypatch.setattr(workflow_real,'build_evidence',lambda *args:ev)
    out=tmp_path/'out'
    monkeypatch.setattr(sys,'argv',['workflow_real','--input-root',str(root),'--output',str(out),
        '--expected-candidates','13','--case','gel','--publication','--cpus','1'])
    workflow_real.main()
    assert len(calls)==1
    evidence=pl.read_parquet(out/'gel/alignment_evidence.parquet')
    assert evidence.equals(ev.alignments)
    provenance=json.loads((out/'provenance.json').read_text())
    parameters=provenance['cases'][0]['error_model_parameters']
    assert parameters['conf_x']==model._conf_x.tolist()
    assert parameters['conf_e']==model._conf_e.tolist()
    assert parameters['substitution']==model.substitution.tolist()
    assert parameters['e_db']==model.e_db
    for aggregation in ['alignment_union','exclusive_spectrum']:
        expected=ev.score(aggregation=aggregation)
        actual=pl.read_parquet(out/'gel'/f'{aggregation}_scores.parquet')
        assert actual.equals(expected)
    layers=pl.read_csv(out/'method_layers.tsv',separator='\t')
    assert layers['method'].to_list()==['match_count','match_count_isotype','no_competition','competition','full']
