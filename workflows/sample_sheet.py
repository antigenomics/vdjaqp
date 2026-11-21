"""Validate and load the observed-analysis sample sheet."""

import csv
import re
from pathlib import Path


def read_samples(path):
    required = {"sample", "candidates", "peptides", "producer", "model"}
    samples = {}
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not required.issubset(reader.fieldnames or []) or len(
            set(reader.fieldnames)
        ) != len(reader.fieldnames):
            raise ValueError("Invalid sample sheet columns")
        for row in reader:
            sample = row["sample"]
            if not re.fullmatch(r"[A-Za-z0-9_-]+", sample or "") or sample in samples:
                raise ValueError(f"Invalid or duplicate sample at line {reader.line_num}")
            if None in row or any(not (row[field] or "").strip() for field in required):
                raise ValueError(f"Missing or extra value at line {reader.line_num}")
            for field in ("candidates", "peptides", "model"):
                source = Path(row[field]).resolve()
                if not source.is_file():
                    raise ValueError(f"Missing {field} file at line {reader.line_num}")
                row[field] = str(source)
            samples[sample] = row
    if not samples:
        raise ValueError("Empty sample sheet")
    return samples
