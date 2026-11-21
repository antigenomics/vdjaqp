"""Loaders for the Ig-Seq shared datasets."""

from __future__ import annotations

from pathlib import Path

import polars as pl

__all__ = [
    "DATASETS",
    "REGION_NAMES",
    "dataset_root",
    "load_clonotypes",
    "load_peptides",
    "load_matches",
    "load_targets",
    "ms_stem_of",
    "ALLELE_COLUMNS",
    "to_gene_level",
    "isotype_of",
]

#: The three mass-spec experiments. All share one repertoire.
DATASETS = ("dataset_single_mabs", "dataset_mixa", "dataset_mixa_gel")

UNANNOTATED_REGION = 9

#: ``region_ids`` values, one per residue of a clonotype ``seq``.
REGION_NAMES = {
    0: "FR1", 1: "CDR1", 2: "FR2", 3: "CDR2",
    4: "FR3", 5: "CDR3", 6: "NDN", 7: "FR4", 8: "C", UNANNOTATED_REGION: "unannotated",
}


def dataset_root(root: str | Path | None = None) -> Path:
    """Directory holding the three dataset folders and ``ground_truth/``."""
    if root is not None:
        return Path(root)
    project = Path(__file__).resolve().parent.parent
    private = project / "private_data" / "igseq_datasets"
    here = private if private.is_dir() else project / "igseq_datasets"
    if not here.is_dir():
        raise FileNotFoundError(
            f"igseq_datasets not found at {here}; pass root= explicitly"
        )
    return here


def _read(dataset: str, name: str, root: str | Path | None) -> pl.DataFrame:
    if dataset not in DATASETS:
        raise ValueError(f"unknown dataset {dataset!r}; expected one of {DATASETS}")
    return pl.read_parquet(dataset_root(root) / dataset / f"{name}.parquet")


def load_clonotypes(dataset: str = "dataset_mixa", root=None) -> pl.DataFrame:
    """The repertoire: 33,052 heavy-chain clonotypes, identical across datasets.

    ``seq`` is the amino-acid sequence V through J plus 13 constant residues.
    ``region_ids`` labels every residue of it (see :data:`REGION_NAMES`).
    """
    return _read(dataset, "clonotypes", root)


def load_peptides(dataset: str = "dataset_mixa", root=None) -> pl.DataFrame:
    """The mass-spec side: one row per LC-MS/MS peptide after quality filtering.

    ``source_file`` names the antibody in ``dataset_single_mabs`` and leaks the
    answer; use :func:`ms_stem_of` to recover which run a peptide came from for
    per-mAb *evaluation*, and never let it reach a scoring path.
    """
    return _read(dataset, "peptides", root)


def load_matches(dataset: str = "dataset_mixa", root=None) -> pl.DataFrame:
    """The *shipped* peptide-to-clonotype join, for reproduction testing only.

    Produced by the upstream seed-and-extend matcher at budget ``floor(L/6)``.
    Used by ``tests/test_match.py`` to prove :mod:`vdjaqp.match` reproduces it,
    and nowhere else.
    """
    return _read(dataset, "matches", root)


def load_targets(root=None) -> pl.DataFrame:
    """The 12 known-answer clonotypes.

    ``is_calibration_reference`` marks the 8 used to fit the de novo error model;
    ``has_single_run`` marks the 8 that also have an individual MS run. Filter
    this table rather than hard-coding 8 or 12.
    """
    return pl.read_csv(dataset_root(root) / "ground_truth" / "target_clonotypes.csv")


def ms_stem_of(peptides: pl.DataFrame) -> pl.Series:
    """Which antibody each ``dataset_single_mabs`` peptide was run from.

    The leading field of ``source_file`` (``reference-TRYP-1.raw`` -> ``reference``). These
    are *ms_stem* values, not clone ids: join through
    :func:`load_targets` to reach ``clone_id``. Evaluation only -- a score that
    sees this is scoring the file name, not the sequence.
    """
    return peptides["source_file"].str.split("-").list.get(0).alias("ms_stem")


#: Metadata columns whose *allele* suffix is an artefact of how the 12 reference
#: antibodies were injected, not biology. See :func:`to_gene_level`.
ALLELE_COLUMNS = ("v_name", "j_name", "c_name")


def to_gene_level(clonotypes: pl.DataFrame) -> pl.DataFrame:
    """Strip allele suffixes from ``v_name``, ``j_name`` and ``c_name``.

    Apply before scoring to prevent injected-reference allele labels from leaking
    ground truth into features or priors. Preserve alternative gene calls in their
    producer order; collapse only alternatives belonging to the same gene. The
    original producer provenance stays with the private input records.
    """
    return clonotypes.with_columns(
        [pl.col(c).str.split(",").list.eval(
            pl.element().str.strip_chars().str.split("*").list.get(0)
         ).list.unique(maintain_order=True).list.join(",").alias(c) for c in ALLELE_COLUMNS
         if c in clonotypes.columns]
    )


def isotype_of(clonotypes: pl.DataFrame) -> pl.Series:
    """IgG / IgM / other / unknown, from the constant-region call.

    Antibodies purified from serum are overwhelmingly IgG while the repertoire
    is 87% IgM, so this is the single most useful prior available. Derived at
    gene level, so it cannot carry the allele leak.
    """
    gene = pl.element().str.strip_chars().str.split("*").list.get(0)
    classes = clonotypes["c_name"].str.split(",").list.eval(
        pl.when(gene.is_in(["IGHG", "IGHG1", "IGHG2", "IGHG3", "IGHG4"])).then(pl.lit("IgG"))
        .when(gene == "IGHM").then(pl.lit("IgM"))
        .when(gene.is_null() | gene.is_in(["", "unknown", "IGHC"])).then(pl.lit("unknown"))
        .otherwise(pl.lit("other"))
    ).list.unique()
    return (
        pl.when(classes.list.len() == 1).then(classes.list.first())
        .otherwise(pl.lit("unknown"))
        .alias("isotype")
    )
