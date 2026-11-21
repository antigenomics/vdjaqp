"""Lossless interchange must preserve observed sequence and anchor coordinates."""
import polars as pl
import pytest

from vdjaqp import adapter


def example():
    return pl.DataFrame({
        "clone_id": ["clone_1"], "seq": ["AACAWWWGG"],
        "cdr3_seq": ["CAWWW"], "cdr3_start": [2], "cdr3_end": [7],
        "region_ids": [[0, 0, 5, 6, 6, 6, 7, 8, 8]],
        "v_name": ["IGHV1-2*01"], "j_name": ["IGHJ4*01"],
        "c_name": ["IGHG1*01"], "entropy_observed": [False],
    })


def test_lossless_observed_sequence_roundtrip():
    assert hasattr(adapter, "to_airr_compatible"), "adapter is not implemented"
    original = example()
    airr = adapter.to_airr_compatible(original, source_format="supplied-export-v1")
    assert airr["sequence_aa"][0] == original["seq"][0]
    assert airr["junction_aa"][0] == "CAWWW"
    assert "cdr3_aa" not in airr.columns
    assert airr["igseq_source_format"][0] == "supplied-export-v1"
    restored = adapter.from_airr_compatible(airr)
    assert restored.select(original.columns).equals(original)


@pytest.mark.parametrize("change", [
    {"region_ids": [[0]]}, {"cdr3_end": [6]}, {"cdr3_seq": [None]},
    {"clone_id": [None]}, {"cdr3_start": [-1]},
    {"region_ids": [[0, 0, 5, None, 6, 6, 7, 8, 8]]},
    {"region_ids": [[0, 0, 5, 10, 6, 6, 7, 8, 8]]},
    {"clone_id": [""]},
])
def test_rejects_incomplete_or_inconsistent_sequence_truth(change):
    frame = example().with_columns(pl.Series(k, v) for k, v in change.items())
    with pytest.raises(ValueError):
        adapter.to_airr_compatible(frame, source_format="supplied-export-v1")


def test_rejects_ambiguous_columns_and_missing_sequence():
    with pytest.raises(ValueError, match="collision"):
        adapter.to_airr_compatible(example().with_columns(sequence_aa=pl.lit("OTHER")), source_format="v1")
    with pytest.raises(ValueError, match="missing"):
        adapter.to_airr_compatible(example().drop("seq"), source_format="v1")
    with pytest.raises(ValueError, match="unique"):
        adapter.to_airr_compatible(pl.concat([example(), example()]), source_format="v1")
    with pytest.raises(ValueError, match="source_format"):
        adapter.to_airr_compatible(example(), source_format="")


def test_reverse_requires_producer_format():
    airr = adapter.to_airr_compatible(example(), source_format="supplied-export-v1")
    with pytest.raises(ValueError, match="source_format"):
        adapter.from_airr_compatible(airr.drop("igseq_source_format"))


def test_cli_converts_without_overwriting(tmp_path):
    import subprocess
    import sys
    source, output = tmp_path / "input.parquet", tmp_path / "adapted.parquet"
    example().write_parquet(source)
    command = [sys.executable, "-m", "vdjaqp.adapter", str(source), str(output),
               "--source-format", "supplied-export-v1"]
    subprocess.run(command, check=True, capture_output=True)
    assert output.exists(), "conversion CLI did not write its output"
    assert pl.read_parquet(output)["junction_aa"][0] == "CAWWW"
    before = output.read_bytes()
    assert subprocess.run(command, capture_output=True).returncode != 0
    assert output.read_bytes() == before
