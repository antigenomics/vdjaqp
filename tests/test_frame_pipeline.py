"""Frame input uses the same search and scoring path as file input."""
import numpy as np
import polars as pl
from vdjaqp.pipeline import Evidence
from vdjaqp.errormodel import ErrorModel, mass_prior


def test_frame_and_file_routes_agree(monkeypatch):
    import vdjaqp.pipeline as pipeline
    import vdjaqp.background as background
    monkeypatch.setattr(background, 'canonical_v', lambda: {})
    monkeypatch.setattr(background, 'canonical_j', lambda: {})
    clones = pl.DataFrame({'clone_id':['c0'], 'seq':['ACDEFGHKLMNP'],
        'region_ids':[[6]*12], 'cdr3_start':[0], 'cdr3_end':[12],
        'v_name':['IGHV1-2'], 'j_name':['IGHJ4'], 'c_name':['IGHG1']})
    peptides = pl.DataFrame({'peptide_id':['run/100/2/50/10/1/ACDEFG'],
        'seq':['ACDEFG'], 'local_confidence':[[float('nan')]*6], 'search_type':['db'],
        'alc':[float('nan')], 'alc_rank':[1.], 'n_reads':[1]})
    error = ErrorModel(np.array([0.,1.]),np.array([.1,.01]),mass_prior(),.1)
    assert hasattr(Evidence,'from_frames'), 'shared frame entry point missing'
    monkeypatch.setattr(pipeline,'load_clonotypes',lambda *a:clones)
    monkeypatch.setattr(pipeline,'load_peptides',lambda *a:peptides)
    monkeypatch.setattr(ErrorModel,'calibrate',lambda **kw:error)
    direct = Evidence.from_frames(clones,peptides,error_model=error,max_subs=0)
    loaded = Evidence.build(max_subs=0)
    assert direct.matches.equals(loaded.matches)
    assert direct.alignments.equals(loaded.alignments)
    assert direct.score().equals(loaded.score())
