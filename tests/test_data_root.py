"""A local checkout finds its ignored inputs without environment variables."""
from vdjaqp import data


def test_private_root_default_and_explicit_override(tmp_path, monkeypatch):
    monkeypatch.setattr(data, "__file__", str(tmp_path / "vdjaqp" / "data.py"))
    private = tmp_path / "private_data" / "igseq_datasets"
    private.mkdir(parents=True)
    assert data.dataset_root() == private
    assert data.dataset_root(tmp_path / "explicit") == tmp_path / "explicit"
