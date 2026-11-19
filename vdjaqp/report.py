"""Per-clonotype immunological readout: what was sequenced, what was measured."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import polars as pl

from .data import REGION_NAMES, isotype_of, load_targets, to_gene_level
from .germline import aligned_columns, canonical_j, canonical_v, common_reference

__all__ = ["ClonotypeProfile", "profile", "clonotype_table", "shm_positions", "peptide_support_table", "placement_depth"]

#: Constant-region label; legacy inputs use an amino-acid stub. Explicit
#: native nucleotide support can also mark observed constant residues.
_CONST = 8

#: ``region_ids`` values inside CDR3: templated flanks and non-templated core.
_JUNCTION = (5, 6)

#: How far past the junction boundaries to align germline when calling somatic
#: hypermutation, matching :mod:`vdjaqp.shm`.
_V_WINDOW, _J_WINDOW = 160, 40


def shm_positions(seq: str, v_gene: str, j_gene: str, regions: np.ndarray,
                  cdr3_start: int, cdr3_end: int) -> np.ndarray:
    """Boolean mask of residues that differ from their germline segment.

    Only V- and J-templated positions can be called: the non-templated core has
    no germline to compare against, so those positions are always ``False``
    rather than silently counted as mutated.
    """
    mutated = np.zeros(len(seq), dtype=bool)
    core = regions == 6

    # aligned_columns yields (query_position, reference_position, reference
    # residue) -- the middle value is an index, not a residue. Compare against
    # the query's own letter, never against that index.
    match = common_reference(v_gene, canonical_v())
    if match is not None:
        _, reference = match
        stop = min(len(seq), cdr3_start + _V_WINDOW)
        for position, _, germline in aligned_columns(seq[:stop], reference):
            if position < len(seq) and not core[position]:
                mutated[position] = seq[position] != germline

    match = common_reference(j_gene, canonical_j())
    if match is not None:
        _, reference = match
        begin = max(0, cdr3_end - _J_WINDOW)
        for position, _, germline in aligned_columns(seq[begin:], reference):
            index = begin + position
            if index < len(seq) and not core[index]:
                mutated[index] = seq[index] != germline
    return mutated


@dataclass(frozen=True)
class ClonotypeProfile:
    """Everything known about one clonotype, aligned to its own residues.

    Every array is one value per residue of :attr:`seq`, so they can be plotted
    directly against each other.
    """

    clone_id: str
    seq: str
    regions: np.ndarray          #: IMGT region id per residue
    entropy: np.ndarray          #: supplied variability, 0 to 1; NaN if unavailable
    nt_support: np.ndarray       #: explicit observed support, or legacy variable-region proxy
    shm: np.ndarray              #: residue differs from its germline segment
    depth: np.ndarray            #: peptides covering the residue
    agreeing_depth: np.ndarray   #: peptides covering it *and* agreeing
    information: np.ndarray      #: nats gained by confirming the residue
    peptides: pl.DataFrame       #: one row per aligning peptide
    is_probe: bool               #: one of the injected reference antibodies

    @property
    def region_names(self) -> list[str]:
        return [REGION_NAMES[int(r)] for r in self.regions]

    def covered(self, regions=_JUNCTION) -> float:
        """Fraction of residues in ``regions`` reached by at least one peptide."""
        mask = np.isin(self.regions, regions)
        return float((self.depth[mask] > 0).mean()) if mask.any() else float("nan")


def profile(clone_id: str, evidence, targets: set[str] | None = None
            ) -> ClonotypeProfile:
    """Build the per-residue profile of one clonotype.

    Args:
        clone_id: Which clonotype.
        evidence: An :class:`vdjaqp.pipeline.Evidence`, already built.
        targets: Injected probe ids; loaded from the ground truth if omitted.
    """
    clonotypes = evidence.clonotypes
    row_index = {c: i for i, c in enumerate(clonotypes["clone_id"])}
    if clone_id not in row_index:
        raise KeyError(f"unknown clone_id {clone_id!r}")
    row = clonotypes.row(row_index[clone_id], named=True)
    gene = to_gene_level(clonotypes).row(row_index[clone_id], named=True)

    seq = row["seq"]
    regions = np.asarray(row["region_ids"])
    entropy = (np.asarray(row["entropy"], dtype=float)
               if row.get("entropy") is not None else np.full(len(seq), np.nan))

    # Explicit native support can include observed constant sequence. The legacy
    # fallback is a supplied-region proxy, not proof of nucleotide ancestry.
    if row.get("nt_support") is None:
        nt_support = regions != _CONST
    else:
        nt_support = np.asarray(row["nt_support"])
        if nt_support.dtype != np.bool_ or nt_support.shape != (len(seq),):
            raise ValueError("nucleotide support must be one Boolean per residue")

    shm = shm_positions(seq, gene["v_name"], gene["j_name"], regions,
                        int(row["cdr3_start"]), int(row["cdr3_end"]))

    hits = (evidence.alignments.filter(pl.col("clone_id") == clone_id)
            .join(evidence.peptides.select(
                      "peptide_id",
                      pl.col("seq").alias("peptide"),
                      "alc", "search_type", "enzyme"),
                  on="peptide_id", how="left")
            .sort("llr", descending=True))

    depth, agreeing = placement_depth(seq, hits)

    offset = int(evidence.background.offset[row_index[clone_id]])
    information = evidence.background.information()[offset:offset + len(seq)]

    probes = targets if targets is not None else set(load_targets()["clone_id"])
    return ClonotypeProfile(
        clone_id=clone_id, seq=seq, regions=regions, entropy=entropy,
        nt_support=nt_support, shm=shm, depth=depth, agreeing_depth=agreeing,
        information=information, peptides=hits, is_probe=clone_id in probes,
    )


def placement_depth(seq: str, hits: pl.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Compatible and agreeing depth from retained placements, with no scoring.

    ``peptide`` contains observation sequences; intervals are zero-based half-open.
    This is the same coverage calculation used by :func:`profile`.
    """
    depth = np.zeros(len(seq), dtype=int)
    agreeing = np.zeros(len(seq), dtype=int)
    for start, end, peptide in zip(hits["start"], hits["end"], hits["peptide"]):
        depth[start:end] += 1
        if peptide is None:
            continue
        window = np.frombuffer(seq[start:end].encode(), dtype=np.uint8)
        read = np.frombuffer(peptide.encode(), dtype=np.uint8)
        n = min(window.size, read.size)
        agreeing[start:start + n] += _same(window[:n], read[:n])
    return depth, agreeing


def _same(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Residue-wise agreement under the MS alphabet, where I and L are one."""
    return ((a == b) | ((a == ord("I")) & (b == ord("L")))
            | ((a == ord("L")) & (b == ord("I")))).astype(int)


def clonotype_table(evidence, clone_ids=None, *, scored: pl.DataFrame | None = None,
                    targets: set[str] | None = None) -> pl.DataFrame:
    """One row per clonotype, in the fields an immunologist reads first.

    Args:
        evidence: A built :class:`vdjaqp.pipeline.Evidence`.
        clone_ids: Which clonotypes; defaults to the injected probes.
        scored: Output of ``evidence.score()``; computed if omitted.
        targets: Probe ids, for the ``is_probe`` column.

    Returns:
        ``clone_id``, gene calls, isotype, junction and IMGT CDR3 with lengths,
        lineage, SHM load, nucleotide length, peptide and coverage counts, and
        the score columns including ``framework_only``.
    """
    probes = targets if targets is not None else set(load_targets()["clone_id"])
    if clone_ids is None:
        clone_ids = [c for c in evidence.clonotypes["clone_id"] if c in probes]
    table = scored if scored is not None else evidence.score()
    ranked = (table.sort(["score", "clone_id"], descending=[True, False])
              .with_row_index("rank", offset=1))
    by_id = {r["clone_id"]: r for r in ranked.iter_rows(named=True)}

    gene = to_gene_level(evidence.clonotypes)
    isotype = dict(zip(gene["clone_id"],
                       evidence.clonotypes.select(isotype_of(evidence.clonotypes))["isotype"]))
    meta = {r["clone_id"]: r for r in gene.iter_rows(named=True)}

    rows = []
    for clone_id in clone_ids:
        p = profile(clone_id, evidence, targets=probes)
        m, s = meta[clone_id], by_id[clone_id]
        junction = m["cdr3_seq"]
        templated = p.regions != 6
        ndn_coverage = p.covered(regions=(6,))
        rows.append({
            "clone_id": clone_id,
            "is_probe": p.is_probe,
            "rank": s["rank"],
            "v_gene": m["v_name"],
            "j_gene": m["j_name"],
            "isotype": isotype[clone_id],
            "lineage_id": m["lineage_id"],
            "junction": junction,
            "cdr3_imgt": junction[1:-1],
            "cdr3_len": max(len(junction) - 2, 0),
            "ndn_len": int((p.regions == 6).sum()),
            "aa_len": len(p.seq),
            "nt_len": len(m["nt_seq"]),
            "shm_aa": int(p.shm.sum()),
            "shm_pct": round(100 * p.shm.sum() / max(templated.sum(), 1), 1),
            "n_peptides": p.peptides.height,
            "ndn_covered_pct": round(100 * ndn_coverage) if np.isfinite(ndn_coverage) else None,
            "confirmed_nats": round(float(s["confirmed_nats"]), 2),
            "junction_nats": round(float(s["junction_nats"]), 2),
            "framework_only": bool(s["framework_only"]),
        })
    return pl.DataFrame(rows).sort("rank")


def peptide_support_table(evidence) -> pl.DataFrame:
    """Count distinct peptide rows/species per candidate from fresh native matches.

    I/L species are collapsed; exact species require zero native substitutions.
    Shared peptides may support several candidates; zero-hit candidates remain.
    """
    support = (evidence.matches.join(
        evidence.peptides.select('peptide_id', pl.col('seq').str.replace_all('I', 'L').alias('species')),
        on='peptide_id', validate='m:1')
        .group_by('clone_id').agg(
            pl.col('peptide_id').n_unique().alias('compatible_peptide_rows'),
            pl.col('species').n_unique().alias('compatible_peptide_species'),
            pl.col('species').filter(pl.col('n_subs') == 0).n_unique().alias('exact_peptide_species')))
    return evidence.clonotypes.select('clone_id').join(support, on='clone_id', how='left', validate='1:1').fill_null(0)
