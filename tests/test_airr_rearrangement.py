import gzip
from pathlib import Path

import pytest

from vdjaqp.airr_rearrangement import read_rearrangements


EXAMPLE = Path(__file__).resolve().parents[1] / "examples/inputs/rearrangements.tsv"


def test_airr_gzip_and_boundaries(tmp_path):
    file = tmp_path / "rearrangements.tsv.gz"
    file.write_bytes(gzip.compress(EXAMPLE.read_bytes()))
    row = read_rearrangements(file).row(0, named=True)
    assert row["seq"] == "ACDEFGHCLMNPWYG"
    assert row["seq"][row["cdr3_start"]:row["cdr3_end"]] == "CLMNPW"
    assert row["region_ids"] == [4] * 8 + [5] * 4 + [7] * 3
    assert row["background_member"] is True
    assert row["c_name"] is None


def test_airr_translates_from_nucleotide_junction_phase(tmp_path):
    file = tmp_path / "rearrangements.tsv"
    file.write_text(EXAMPLE.read_text().replace("ACDEFGHC.LMNPWYG", "")
                    .replace("\tCLMNPW\t", "\t\t").replace("GCTTGT", "GCT...TGT"))
    assert read_rearrangements(file)["seq"].item() == "ACDEFGHCLMNPWYG"


@pytest.mark.parametrize(("old", "new", "message"), [
    ("\tLMNP\t", "\tCLMNPW\t", "without its anchors"),
    ("ACDEFGHC.LMNPWYG", "ACDEFGHC.LMNPWCLMNPWYG", "exactly once"),
])
def test_airr_rejects_inconsistent_sequence(tmp_path, old, new, message):
    file = tmp_path / "rearrangements.tsv"
    file.write_text(EXAMPLE.read_text().replace(old, new))
    with pytest.raises(ValueError, match=message):
        read_rearrangements(file)


def test_airr_preserves_unknown_regions_and_explicit_roles(tmp_path):
    header, row = EXAMPLE.read_text().splitlines()
    file = tmp_path / "rearrangements.tsv"
    file.write_text(header + "\tbackground_member\n"
                    + row.replace("\tACDEFGHC\t", "\t\t") + "\tF\n")
    imported = read_rearrangements(file).row(0, named=True)
    assert imported["region_ids"] == [9] * 8 + [5] * 4 + [7] * 3
    assert imported["background_member"] is False
