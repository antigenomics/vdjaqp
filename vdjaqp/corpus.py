"""Known antibodies, and whether a donor clonotype is one of them."""

from __future__ import annotations

import ast
from pathlib import Path

import polars as pl

__all__ = ["CORPUS_COLUMNS", "load_sabdab", "load_plabdab", "load_corpus",
           "match_corpus"]

#: One row per known antibody. ``cdr3_imgt`` is the join key: IMGT CDR3, anchors
#: excluded, which is what both corpora publish and what a donor junction
#: becomes when its first and last residue are dropped.
CORPUS_COLUMNS = ("source", "antibody_id", "cdr3_imgt", "antigen",
                  "antigen_species", "vh")


def _external(root: str | Path | None) -> Path:
    if root is not None:
        return Path(root)
    here = Path(__file__).resolve().parent.parent / "external"
    if not here.is_dir():
        raise FileNotFoundError(f"external corpora not found at {here}")
    return here


def load_sabdab(root=None) -> pl.DataFrame:
    """SAbDab: structural antibodies with an antigen name and species.

    ``CDR-H3`` is IMGT CDR3 -- anchors excluded, median length 14 against the
    donor junction's 16.
    """
    frame = pl.read_csv(
        _external(root) / "sabdab_summary_all.csv",
        infer_schema_length=5000, ignore_errors=True,
    )
    return (
        frame.select(
            pl.lit("sabdab").alias("source"),
            pl.col("PDB").cast(pl.String).alias("antibody_id"),
            pl.col("CDR-H3").alias("cdr3_imgt"),
            pl.col("antigen_name").alias("antigen"),
            pl.col("antigen_species").cast(pl.String).alias("antigen_species"),
            pl.col("VH").alias("vh"),
        )
        .drop_nulls("cdr3_imgt")
        .filter(pl.col("cdr3_imgt").str.len_chars() > 0)
        .unique()
    )


def load_plabdab(root=None) -> pl.DataFrame:
    """PLAbDab's nanobody (VHH) set. CDRs arrive as a stringified dict."""
    frame = pl.read_csv(_external(root) / "plabdab_nano" / "vhh_sequences.csv.gz")

    def _cdrh3(cell: str | None) -> str | None:
        if not cell:
            return None
        try:
            return ast.literal_eval(cell).get("CDRH3")
        except (ValueError, SyntaxError):
            return None

    return (
        frame.with_columns(
            pl.col("cdr_sequences")
            .map_elements(_cdrh3, return_dtype=pl.String)
            .alias("cdr3_imgt")
        )
        .select(
            pl.lit("plabdab_nano").alias("source"),
            pl.col("ID").cast(pl.String).alias("antibody_id"),
            pl.col("cdr3_imgt"),
            pl.col("targets_mentioned").cast(pl.String).alias("antigen"),
            pl.col("organism").cast(pl.String).alias("antigen_species"),
            pl.col("sequence").alias("vh"),
        )
        .drop_nulls("cdr3_imgt")
        .unique()
    )


def load_corpus(root=None) -> pl.DataFrame:
    """Every known antibody available locally, as one table of :data:`CORPUS_COLUMNS`.

    QuantAb (37,359 unique heavy chains with antigen and Kd) and AbDb (78,573
    structures) are not included: neither publishes a CDR3, so both need an arda
    annotation pass first.
    """
    return pl.concat(
        [load_sabdab(root).select(CORPUS_COLUMNS),
         load_plabdab(root).select(CORPUS_COLUMNS)],
        how="vertical",
    )


def match_corpus(clonotypes: pl.DataFrame, corpus: pl.DataFrame | None = None,
                 root=None) -> pl.DataFrame:
    """Match donor clonotypes into the corpus on exact IMGT CDR3.

    Args:
        clonotypes: Must carry ``clone_id`` and ``cdr3_seq``, the latter being
            the **junction** (anchors included), as the shipped datasets provide it.
        corpus: Output of :func:`load_corpus`; loaded if omitted.
        root: ``external/`` location.

    Returns:
        One row per (clonotype, corpus entry) match with ``clone_id``,
        ``junction``, ``cdr3_imgt``, ``source``, ``antibody_id``, ``antigen``,
        ``antigen_species``. Join it to :func:`vdjaqp.pgen.junction_pgen` before
        believing any of it.
    """
    corpus = corpus if corpus is not None else load_corpus(root)
    donor = clonotypes.select(
        "clone_id",
        pl.col("cdr3_seq").alias("junction"),
        pl.col("cdr3_seq")
        .str.slice(1, pl.col("cdr3_seq").str.len_chars() - 2)
        .alias("cdr3_imgt"),
    ).filter(pl.col("cdr3_imgt").str.len_chars() > 0)
    return donor.join(corpus, on="cdr3_imgt", how="inner")
