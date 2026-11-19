"""The background model: what residue would we expect at this position anyway?"""

from __future__ import annotations

import numpy as np
import polars as pl
from .alphabet import N_AA, encode_frame
from .data import to_gene_level
from .germline import aligned_columns, canonical_j, canonical_v, common_reference

__all__ = ["Background"]

#: Dirichlet pseudocount pulling a cell towards its backoff level.
_LAMBDA = 20.0

#: How far past each CDR3 boundary a germline alignment may reach, so the V-
#: and J-templated flanks of CDR3 are anchored rather than pooled into the
#: non-templated core.
_V_OVERHANG, _J_OVERHANG = 4, 6

#: ``region_ids`` value for the non-templated junction core.
_NDN = 6

#: ``region_ids`` value for the appended constant-region stub. Those residues
#: are a fixed-length copy of the isotype's CH1 start, so position within the
#: block is an exact coordinate -- they must not fall through to the junction
#: cell, which would price an invariant constant residue at 2.1 nats and
#: over-credit every peptide spanning the J-C junction.
_CONSTANT = 8


class Background:
    """Leave-one-out residue frequencies at germline-anchored coordinates.

    Everything is held as per-*cell* tables of shape ``(n_cells, 19)`` plus a
    flat ``cell_id`` per residue. The leave-one-out correction touches only the
    clonotype's own letter, so a dense ``(4.4M, 19)`` frequency array is never
    materialised -- the tables below are a few hundred kilobytes.

    Attributes:
        cell_id: Cell of every residue, concatenated over clonotypes.
        code: Alphabet index of every residue, ``-1`` outside the alphabet.
        offset: Where each clonotype's residues begin.
        length: Residues per clonotype.

    Optional nonmissing Boolean ``background_member`` supplies the source
    population role independently of candidate identifiers. Without that column,
    the legacy ``exclude_prefix`` rule is retained for compatibility.
    """

    def __init__(self, clonotypes: pl.DataFrame, exclude_prefix: str = "Ab_"):
        clonotypes = to_gene_level(clonotypes)
        self.code, self.offset, self.length = encode_frame(clonotypes["seq"])
        self.cell_id, self.n_cells = self._anchor(clonotypes)
        if not hasattr(self, "context_valid"):
            self.context_valid = np.ones(self.code.size, dtype=bool)
        self._tabulate(clonotypes, exclude_prefix)

    # -- coordinates -----------------------------------------------------
    def _anchor(self, clonotypes: pl.DataFrame) -> tuple[np.ndarray, int]:
        """Give every residue the germline position it aligns to.

        Residues aligning to a gap are insertions relative to germline, hence
        non-templated, and fall through to the junction cell -- as does
        everything labelled NDN. Missing or unresolved gene context is recorded
        separately as U and supplies no positional identification evidence.
        """
        v_germ, j_germ = canonical_v(), canonical_j()
        # Cell 0 is the non-templated junction; every other cell is one
        # (side, gene, germline position) triple, recorded for the backoff.
        cell_index = {("N", "", -1): 0}

        cell_id = np.zeros(self.code.size, dtype=np.int32)
        seqs = clonotypes["seq"].to_list()
        regions = clonotypes["region_ids"].to_list()
        starts = clonotypes["cdr3_start"].to_numpy()
        ends = clonotypes["cdr3_end"].to_numpy()
        v_genes = clonotypes["v_name"].to_list()
        j_genes = clonotypes["j_name"].to_list()
        c_genes = clonotypes["c_name"].to_list()
        self.n_unanchored = 0

        for row, seq in enumerate(seqs):
            base = self.offset[row]
            block = cell_id[base:base + self.length[row]]
            v_match = common_reference(v_genes[row], v_germ)
            if v_match is None:
                self.n_unanchored += 1
            else:
                v_gene, v_ref = v_match
                stop = min(len(seq), int(starts[row]) + _V_OVERHANG)
                self._write_positions(block, seq[:stop], v_ref, 0, "V",
                                      v_gene, cell_index)
            j_match = common_reference(j_genes[row], j_germ)
            if j_match is not None:
                j_gene, j_ref = j_match
                begin = max(0, int(ends[row]) - _J_OVERHANG)
                self._write_positions(block, seq[begin:], j_ref, begin, "J",
                                      j_gene, cell_index)
            region = np.asarray(regions[row])
            unresolved = ((region <= 4) & (v_match is None)) | (
                (region == 7) & (j_match is None))
            if v_match is None or j_match is None:
                unresolved |= (region == 5) & (block == 0)
            if unresolved.any():
                key = ("U", "", -1)
                if key not in cell_index:
                    cell_index[key] = len(cell_index)
                block[unresolved] = cell_index[key]
            block[region == _NDN] = 0
            constant = np.flatnonzero(region == _CONSTANT)
            for offset, at in enumerate(constant):
                key = ("C", c_genes[row], offset)
                known = cell_index.get(key)
                if known is None:
                    known = cell_index[key] = len(cell_index)
                block[at] = known

        self.cells = pl.DataFrame(
            {
                "cell_id": np.fromiter(cell_index.values(), dtype=np.int32,
                                       count=len(cell_index)),
                "side": [k[0] for k in cell_index],
                "gene": [k[1] for k in cell_index],
                "germline_pos": [k[2] for k in cell_index],
            },
            schema={"cell_id": pl.Int32, "side": pl.String, "gene": pl.String,
                    "germline_pos": pl.Int32},
        ).sort("cell_id")
        self.context_valid = np.asarray(self.cells["side"] != "U")[cell_id]
        return cell_id, len(cell_index)

    @staticmethod
    def _write_positions(block, query, reference, shift, side, gene, cell_index):
        """Record germline positions for one aligned segment."""
        for query_at, reference_at, _ in aligned_columns(query, reference):
            key = (side, gene, reference_at)
            known = cell_index.get(key)
            if known is None:
                known = cell_index[key] = len(cell_index)
            block[shift + query_at] = known

    # -- tables ----------------------------------------------------------
    def _tabulate(self, clonotypes: pl.DataFrame, exclude_prefix: str) -> None:
        """Count residues per cell, then build the smoothed numerator."""
        if "background_member" in clonotypes.columns:
            membership = clonotypes["background_member"]
            if membership.dtype != pl.Boolean or membership.null_count():
                raise ValueError("background_member must be Boolean and nonmissing")
            excluded = ~membership.to_numpy()
        else:
            excluded = clonotypes["clone_id"].str.starts_with(exclude_prefix).to_numpy()
        keep = ~np.repeat(excluded, self.length) & (self.code >= 0)
        pooled = np.bincount(self.code[keep], minlength=N_AA).astype(float)
        keep &= getattr(self, "context_valid", True)
        self._contributed = keep

        counts = np.zeros((self.n_cells, N_AA), dtype=np.float64)
        np.add.at(counts, (self.cell_id[keep], self.code[keep]), 1.0)

        if not pooled.sum():
            raise ValueError("background requires at least one canonical residue from a background member")
        pooled /= pooled.sum()
        backoff = self._backoff(counts, pooled)

        #: Smoothed count numerator; leave-one-out subtracts 1 from the
        #: clonotype's own letter at query time.
        self.numerator = counts + _LAMBDA * backoff
        self.total = counts.sum(axis=1)
        self.pooled = pooled

    def _backoff(self, counts: np.ndarray, pooled: np.ndarray) -> np.ndarray:
        """Per-cell prior: V gene backs off to V family, J and junction to pooled.

        A cell holding only a handful of clonotypes -- a rare V gene, a germline
        position past the end of most sequences -- would otherwise report its
        few observations as certainty.
        """
        family = (
            self.cells.with_columns(
                _families=pl.col("gene").str.split(",").list.eval(
                    pl.element().str.split("-").list.get(0)).list.unique()
            ).with_columns(
                group=pl.when((pl.col("side") == "V") & (pl.col("_families").list.len() == 1))
                .then(
                    pl.concat_str(
                        [pl.col("_families").list.first(),
                         pl.col("germline_pos").cast(pl.String)],
                        separator=":",
                    )
                )
                .otherwise(pl.lit("pooled"))
            )
            .with_columns(group_id=pl.col("group").rank("dense").cast(pl.Int32) - 1)
            .sort("cell_id")["group_id"]
            .to_numpy()
        )
        n_groups = int(family.max()) + 1
        group_counts = np.zeros((n_groups, N_AA), dtype=np.float64)
        np.add.at(group_counts, family, counts)
        group_freq = (group_counts + _LAMBDA * pooled) / (
            group_counts.sum(axis=1, keepdims=True) + _LAMBDA
        )
        return group_freq[family]

    # -- query -----------------------------------------------------------
    def frequency(self, residue: np.ndarray, letter: np.ndarray) -> np.ndarray:
        """Leave-one-out ``P(letter)`` at each of the given residue slots.

        Args:
            residue: Indices into the flat residue arrays.
            letter: Alphabet index whose probability is wanted, same shape.
        """
        cell = self.cell_id[residue]
        own = self.code[residue]
        contributed = self._contributed[residue]
        numerator = self.numerator[cell, letter] - ((own == letter) & contributed)
        return np.maximum(numerator, 1e-9) / (self.total[cell] - contributed + _LAMBDA)

    def mixed(self, residue: np.ndarray, letter: np.ndarray,
              substitution: np.ndarray) -> np.ndarray:
        """``sum_b P(b) * S(b -> letter)``: the background seen through an error."""
        cell = self.cell_id[residue]
        own = self.code[residue]
        pooled_mix = self.numerator @ substitution
        contributed = self._contributed[residue]
        numerator = pooled_mix[cell, letter] - substitution[own, letter] * contributed
        return np.maximum(numerator, 1e-9) / (self.total[cell] - contributed + _LAMBDA)

    def information(self) -> np.ndarray:
        """Surprisal, in nats, of each residue's own letter under the background.

        This is the value of confirming that residue: near zero in framework,
        several nats in the non-templated junction.
        """
        slots = np.arange(self.code.size)
        valid = (self.code >= 0) & self.context_valid
        out = np.zeros(self.code.size)
        out[valid] = -np.log(
            self.frequency(slots[valid], self.code[valid])
        )
        return out
