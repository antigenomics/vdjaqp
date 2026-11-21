"""Lossless AIRR-style interchange of an already annotated sequence export.

This is not an AIRR validator or a native ARDA annotator. ``sequence_aa`` must
contain the supplied full sequence; junction-only clone tables are insufficient.
The adapter preserves data, including reconstructed/stub sequence when present; it
does not certify experimentally measured origin. Original
amino-acid coordinates are zero-based, half-open and retained in extension
fields, never presented as AIRR nucleotide coordinates. Extra columns survive.
"""
from __future__ import annotations

import polars as pl

from .data import REGION_NAMES

__all__ = ["to_airr_compatible", "from_airr_compatible"]

_NAMES = {
    "clone_id": "sequence_id", "seq": "sequence_aa", "nt_seq": "sequence",
    "cdr3_seq": "junction_aa", "v_name": "v_call", "j_name": "j_call",
    "c_name": "c_call", "cdr3_start": "igseq_junction_start_aa",
    "cdr3_end": "igseq_junction_end_aa", "region_ids": "igseq_region_ids",
}
_REQUIRED = set(_NAMES) - {"nt_seq"}


def _validate(frame: pl.DataFrame) -> None:
    missing = _REQUIRED - set(frame.columns)
    if missing:
        raise ValueError(f"missing supplied sequence/annotation fields: {sorted(missing)}")
    if frame["clone_id"].n_unique() != frame.height or frame["clone_id"].null_count():
        raise ValueError("clone_id must be non-null and unique")
    for name in ("clone_id", "seq", "cdr3_seq", "v_name", "j_name", "c_name"):
        if frame.schema[name] != pl.String:
            raise ValueError(f"{name} must be a string column")
    for name in ("cdr3_start", "cdr3_end"):
        if not frame.schema[name].is_integer():
            raise ValueError(f"{name} must be an integer column")
    dtype = frame.schema["region_ids"]
    if dtype.base_type() != pl.List or not dtype.inner.is_integer():
        raise ValueError("region_ids must be a list of integers")
    bad = frame.select(pl.any_horizontal(
        pl.col("clone_id").str.len_chars() == 0,
        pl.col("seq").str.len_chars() == 0,
        pl.col("cdr3_start") < 0,
        pl.col("cdr3_end") <= pl.col("cdr3_start"),
        pl.col("cdr3_end") > pl.col("seq").str.len_chars(),
        pl.col("region_ids").list.len() != pl.col("seq").str.len_chars(),
        pl.col("region_ids").list.eval(
            pl.element().is_null() | ~pl.element().is_in(list(REGION_NAMES))
        ).list.any(),
        pl.col("seq").str.slice(pl.col("cdr3_start"), pl.col("cdr3_end") - pl.col("cdr3_start")) != pl.col("cdr3_seq"),
        pl.col("cdr3_seq").is_null(),
    ).fill_null(True).any()).item()
    if bad:
        raise ValueError("inconsistent supplied sequence, junction coordinates or region coverage")


def to_airr_compatible(frame: pl.DataFrame, *, source_format: str) -> pl.DataFrame:
    """Rename export fields without reannotation or changing producer provenance."""
    if not isinstance(source_format, str) or not source_format.strip():
        raise ValueError("source_format must identify the actual input format")
    _validate(frame)
    mapping = {k: v for k, v in _NAMES.items() if k in frame.columns}
    if set(mapping.values()) & set(frame.columns) or "igseq_source_format" in frame.columns:
        raise ValueError("column collision in AIRR-compatible conversion")
    return frame.rename(mapping).with_columns(igseq_source_format=pl.lit(source_format))


def from_airr_compatible(frame: pl.DataFrame) -> pl.DataFrame:
    """Restore scorer columns from interchange retaining supplied sequence coverage.

    Native ARDA junction-only output must first be joined to its observed full
    sequences/annotations; this function never fills them from germline.
    """
    if "igseq_source_format" not in frame.columns or frame.schema["igseq_source_format"] != pl.String:
        raise ValueError("source_format provenance must accompany the sequence export")
    if frame.select((pl.col("igseq_source_format").str.strip_chars() == "").fill_null(True).any()).item():
        raise ValueError("source_format provenance must not be empty")
    mapping = {v: k for k, v in _NAMES.items() if v in frame.columns}
    if set(mapping.values()) & set(frame.columns):
        raise ValueError("column collision in scorer conversion")
    restored = frame.rename(mapping)
    _validate(restored)
    return restored


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="supplied annotated repertoire Parquet")
    parser.add_argument("output", help="new AIRR-style Parquet; must not exist")
    parser.add_argument("--source-format", required=True, help="actual producer/export format")
    args = parser.parse_args()
    converted = to_airr_compatible(pl.read_parquet(args.input), source_format=args.source_format)
    with open(args.output, "xb") as stream:
        converted.write_parquet(stream, compression="zstd")
