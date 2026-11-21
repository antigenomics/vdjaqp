import csv
import gzip
from pathlib import Path


FIXTURES = Path(__file__).resolve().parents[1] / "examples" / "inputs"


def test_mock_ms_fixture():
    with gzip.open(FIXTURES / "mock-ms.tsv.gz", "rt") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len({row["peptide_id"] for row in rows}) == len(rows) == 10_000


def test_tagged_umi_fastq_fixture():
    path = FIXTURES / "sc5p_v2_hs_PBMC_1k_t_cells1pct.fq.gz"
    with gzip.open(path, "rt") as stream:
        records = list(zip(stream, stream, stream, stream, strict=True))
    assert len(records) == 18_036
    for name, seq, plus, qual in records:
        assert name.startswith("@") and plus.startswith("+")
        assert "RX:Z:" in name and "CB:Z:" in name
        assert len(seq.rstrip()) == len(qual.rstrip())
