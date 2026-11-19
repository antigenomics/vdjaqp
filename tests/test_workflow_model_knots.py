"""Keep the external-model loader compatible with native isotonic knots."""
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import polars as pl
import pytest


def test_loader_preserves_native_tied_knots(tmp_path, monkeypatch):
    from vdjaqp.errormodel import mass_prior
    from vdjaqp.pipeline import Evidence
    import vdjaqp.adapter
    import vdjaqp.data

    path = Path(__file__).resolve().parents[1] / 'scripts/workflow_analysis.py'
    spec = importlib.util.spec_from_file_location('workflow_model_test', path)
    workflow = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(workflow)
    model = dict(conf_x=[.2, .5, .5, .9], conf_e=[.5, .4, .2, .1],
                 substitution=mass_prior().tolist(), e_db=.04)
    model_path = tmp_path / 'model.json'
    model_path.write_text(json.dumps(model))
    clones = pl.DataFrame({'clone_id': ['observed'], 'background_member': [True]})
    monkeypatch.setattr(pl, 'read_parquet', lambda _: clones)
    monkeypatch.setattr(vdjaqp.adapter, 'to_airr_compatible', lambda frame, **_: frame)
    monkeypatch.setattr(vdjaqp.adapter, 'from_airr_compatible', lambda frame: frame)
    monkeypatch.setattr(vdjaqp.data, 'to_gene_level', lambda frame: frame)
    monkeypatch.setattr(sys, 'argv', [str(path), str(tmp_path/'candidates'),
                        str(tmp_path/'peptides'), str(model_path),
                        str(tmp_path/'output'), '--producer', 'arda', '--cpus', '1'])

    class ModelReachedEvidence(Exception):
        pass

    def inspect_model(clones, peptides, *, error_model):
        query = np.array([.1, np.nextafter(.5, 0), .5,
                          np.nextafter(.5, 1), 1.])
        np.testing.assert_array_equal(error_model._conf_x, model['conf_x'])
        np.testing.assert_array_equal(error_model._conf_e, model['conf_e'])
        np.testing.assert_array_equal(error_model.error_rate(query),
                                      np.interp(query, model['conf_x'], model['conf_e']))
        raise ModelReachedEvidence

    monkeypatch.setattr(Evidence, 'from_frames', staticmethod(inspect_model))
    with pytest.raises(ModelReachedEvidence):
        workflow.main()

    for invalid in ([.2, .6, .5, .9], [.2, float('nan'), .5, .9]):
        model['conf_x'] = invalid
        model_path.write_text(json.dumps(model))
        with pytest.raises(SystemExit, match='2'):
            workflow.main()
