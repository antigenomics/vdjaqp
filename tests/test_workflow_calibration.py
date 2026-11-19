"""Shared calibration routing; no scientific fitting or peptide search."""
from pathlib import Path
from scripts import workflow_real
from vdjaqp.pipeline import Evidence
from vdjaqp.errormodel import ErrorModel
from vdjaqp import data


def test_common_calibration_uses_joint_root_and_own_candidates(monkeypatch):
    model=object();clones=object();peptides=object();calls={}
    monkeypatch.setattr(ErrorModel,'calibrate',lambda **kw:calls.update(calibration=kw) or model)
    monkeypatch.setattr(data,'load_clonotypes',lambda dataset,root:calls.update(candidates=(dataset,root)) or clones)
    monkeypatch.setattr(Evidence,'from_frames',lambda c,p,**kw:calls.update(frames=(c,p),options=kw) or 'built')
    result=workflow_real.build_evidence('gel',Path('own'),Path('joint'),peptides,'heldout')
    assert result=='built'
    assert calls['calibration']==dict(root=Path('joint'),exclude={'heldout'})
    assert calls['candidates']==('gel',Path('own'))
    assert calls['frames']==(clones,peptides) and calls['options']['error_model'] is model


def test_default_calibration_route_is_unchanged(monkeypatch):
    calls={}
    monkeypatch.setattr(Evidence,'build',lambda dataset,**kw:calls.update(dataset=dataset,**kw) or 'original')
    assert workflow_real.build_evidence('gel',Path('own'),None,None,None)=='original'
    assert calls==dict(dataset='gel',root=Path('own'),peptides=None,exclude=())
