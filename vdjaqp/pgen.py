"""Generation probability of a junction, as its information content."""

from __future__ import annotations

import numpy as np
import polars as pl

from .resources import bounded_threads

__all__ = ["DEFAULT_SOURCE", "junction_pgen", "insertion_pgen",
           "sharing_probability"]

#: Bundled IGH model to use. ``olga`` is the published germline-templated model;
#: ``learned`` and ``arda`` add tandem D at the cost of a different fit.
DEFAULT_SOURCE = "olga"


def junction_pgen(junctions, *, locus: str = "IGH", source: str = DEFAULT_SOURCE,
                  collapse: bool = True, threads: int | None = None) -> pl.DataFrame:
    """Generation probability of each junction, in one batched call.

    Args:
        junctions: Amino-acid junctions, anchors **included** (``C`` ... ``W``).
            Passing IMGT CDR3 instead silently returns zeros.
        locus: Bundled model locus; ``IGH`` for antibody heavy chains.
        source: Bundled model source, ``olga`` / ``learned`` / ``arda``.
        collapse: Collapse each gene to one representative allele. Keep this on:
            60 of 62 IGH V alleles are affected by truncated deletion grids in
            the published model and some, such as ``IGHV4-30-4*01``, have
            identically zero Pgen.
        threads: Worker threads; defaults to this process's CPU budget. The call
            releases the GIL and threads internally, so it is never wrapped in a
            pool of its own.

    Returns:
        ``junction``, ``pgen``, ``neglog10_pgen`` (null where Pgen is 0) and
        ``scoreable``. Row order matches the input.
    """
    threads = bounded_threads(threads)
    from vdjtools.model import native
    from vdjtools.model.bundled import load_bundled

    sequences = list(junctions)
    model = load_bundled(locus, source=source, collapse=collapse)
    values = np.asarray(
        native.pgen_aa_batch(model, sequences,
                             threads=threads),
        dtype=np.float64,
    )
    scoreable = values > 0.0
    neglog = np.where(scoreable, -np.log10(np.where(scoreable, values, 1.0)), np.nan)
    return pl.DataFrame(
        {
            "junction": sequences,
            "pgen": values,
            "neglog10_pgen": neglog,
            "scoreable": scoreable,
        }
    ).with_columns(
        pl.when(pl.col("scoreable")).then(pl.col("neglog10_pgen"))
        .otherwise(None).alias("neglog10_pgen")
    )


def sharing_probability(pgen, n_rearrangements: float) -> np.ndarray:
    """Chance of seeing a junction of this Pgen again in ``n`` independent draws.

    ``1 - (1 - p)^n``. This is what turns a sequence match into a claim: a
    junction found in two donors is evidence they responded to the same antigen
    only when this probability is small. For a repertoire of ``n = 1e6``
    rearrangements, ``Pgen = 1e-7`` recurs with probability 0.095 -- unremarkable
    -- while ``Pgen = 1e-20`` does not recur at all.

    Args:
        pgen: Generation probabilities.
        n_rearrangements: Number of independent rearrangements sampled. Use the
            repertoire's clonotype count, or a cohort's total, not a read count.

    Returns:
        Probability of at least one independent recurrence, same shape as ``pgen``.
    """
    values = np.asarray(pgen, dtype=np.float64)
    depth = np.asarray(n_rearrangements, dtype=np.float64)
    if (depth.ndim != 0 or not np.isfinite(depth) or depth < 0
            or not np.isfinite(values).all() or ((values < 0) | (values > 1)).any()):
        raise ValueError("probabilities must be finite in [0,1] and depth finite and nonnegative")
    if depth == 0:
        return np.zeros_like(values)
    # Certainty is a valid endpoint; log(0) represents its exact limiting value.
    with np.errstate(divide="ignore"):
        return -np.expm1(depth * np.log1p(-values))


#: Bases an N addition can take. Terminal deoxynucleotidyl transferase adds
#: them without a template, so each is one of four at roughly equal odds.
_BASES = 4.0


def insertion_pgen(clonotypes: pl.DataFrame, *, species: str = "human",
                   threads: int | None = None) -> pl.DataFrame:
    """A model-free generation probability from N-addition length alone.

    The non-templated bases are the only part of a rearrangement drawn at
    random: V, D and J contribute germline, and exonuclease trimming picks from
    a handful of lengths, but each N base is one of four. So

        Pgen ~ 4 ** -(len(np1) + len(np2))

    with the D segment excluded, because D is templated and its bases are a
    gene choice rather than a free draw. This is a floor on improbability, not
    a probability: it ignores the V, D and J choices and the trimming, so it
    is systematically higher than the full model's Pgen and is useful as an
    *ordering*, not as a number to quote.

    Its virtue is that it barely moves under somatic hypermutation. The full
    DP requires the templated flanks to match germline exactly, so a mutation
    there depresses or zeroes Pgen; the N-addition *lengths* are set by where
    the V and J boundaries fall and survive a few substitutions either side.
    That makes this the robust estimator for mutated B-cell junctions, and the
    cross-check on :func:`junction_pgen`.

    Args:
        clonotypes: Needs ``clone_id``, ``cdr3_seq`` (the junction), ``v_name``
            and ``j_name``. Gene-level calls are fine -- the annotator resolves
            alleles itself.
        species: Passed through to the junction annotator.
        threads: Defaults to this process's CPU budget.

    Returns:
        ``clone_id``, ``n_insert`` (non-templated bases), ``d_call``,
        ``insertion_pgen`` and ``neglog10_insertion_pgen``.
    """
    threads = bounded_threads(threads)
    from vdjtools.model.junction import annotate_junctions

    from .data import to_gene_level

    frame = to_gene_level(clonotypes)
    annotated = annotate_junctions(
        frame["cdr3_seq"].to_list(),
        frame["v_name"].to_list(),
        frame["j_name"].to_list(),
        species=species,
        threads=threads,
    )
    n_insert = (
        annotated["np1"].fill_null("").str.len_chars()
        + annotated["np2"].fill_null("").str.len_chars()
    )
    return frame.select("clone_id").with_columns(
        n_insert=n_insert.cast(pl.Int32),
        d_call=annotated["d_call"],
        insertion_pgen=np.power(_BASES, -n_insert.to_numpy().astype(np.float64)),
        neglog10_insertion_pgen=n_insert.to_numpy().astype(np.float64)
        * np.log10(_BASES),
    )
