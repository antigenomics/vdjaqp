import csv

import pytest

from workflows.sample_sheet import read_samples


@pytest.mark.parametrize("sample_names", [[], ["a", "a"], ["../a"], ["sample_1"]])
def test_sample_sheet(tmp_path, monkeypatch, sample_names):
    monkeypatch.chdir(tmp_path)
    for name in ("candidates, raw.parquet", "peptides.parquet", "model.json"):
        (tmp_path / name).touch()
    sheet = tmp_path / "samples.csv"
    with sheet.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["sample", "candidates", "peptides", "producer", "model"])
        for sample in sample_names:
            writer.writerow([sample, "candidates, raw.parquet", "peptides.parquet",
                             "producer's export", "model.json"])
    if sample_names == ["sample_1"]:
        row = read_samples(sheet)["sample_1"]
        assert row["candidates"] == str(tmp_path / "candidates, raw.parquet")
        assert row["producer"] == "producer's export"
        (tmp_path / "model.json").unlink()
        with pytest.raises(ValueError, match="Missing model"):
            read_samples(sheet)
    else:
        with pytest.raises(ValueError, match="Empty|Invalid or duplicate"):
            read_samples(sheet)
