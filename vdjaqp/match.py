"""Peptide-to-clonotype search, on top of :class:`seqtree.TextIndex`."""

from __future__ import annotations

import numpy as np
import polars as pl
import pyarrow as pa
import seqtree

from .resources import bounded_threads

__all__ = ["MS_ALPHABET_MAP", "MIN_QUERY_LEN", "ms_collapse", "PeptideMatcher",
           "shipped_budget"]

#: Residues mass spectrometry cannot tell apart, collapsed to one letter before
#: indexing. Leucine and isoleucine are exactly isobaric (113.08406 Da), so a
#: "mismatch" between them is an artefact of the alphabet, not an observation.
#: Merely *near*-isobaric pairs are deliberately left as substitutions and
#: priced by the error model instead.
MS_ALPHABET_MAP = {"I": "L"}

#: Shortest query ``TextIndex`` can answer at ``k=4``; shorter peptides raise.
MIN_QUERY_LEN = 6


def ms_collapse(seqs) -> list[str]:
    """Project sequences into the mass-spectrometry alphabet (I -> L)."""
    return [s.replace("I", "L") for s in seqs]


def shipped_budget(length: int) -> int:
    """The upstream matcher's substitution budget, ``floor(L / 6)``.

    Depends on peptide length alone -- it ignores the per-residue de novo
    confidence that actually predicts error. Kept only to reproduce the shipped
    join in ``tests/test_match.py``.
    """
    return length // 6


class PeptideMatcher:
    """Substring search of peptides into clonotypes.

    Args:
        clone_ids: Clonotype ids, positionally aligned with ``seqs``.
        seqs: Clonotype amino-acid sequences, *not* yet collapsed.
        k: Seed width. ``k`` belongs to the index, not the query, so one build
            answers every length; performance crosses over at ``L = 2k``.

    The index is built over the I->L collapsed text, so queries must be
    collapsed too -- :meth:`search` does that for you.
    """

    def __init__(self, clone_ids, seqs, k: int = 4):
        self.clone_ids = list(clone_ids)
        collapsed = ms_collapse(seqs)
        if len(self.clone_ids) != len(collapsed):
            raise ValueError("clone IDs and sequences must have the same length")
        self._index = seqtree.TextIndex.build(collapsed, "aa", k=k)
        self.k = k

    @classmethod
    def from_clonotypes(cls, clonotypes: pl.DataFrame, k: int = 4) -> "PeptideMatcher":
        """Build from a :func:`vdjaqp.data.load_clonotypes` frame."""
        return cls(clonotypes["clone_id"].to_list(), clonotypes["seq"].to_list(), k=k)

    def search(
        self,
        peptide_ids,
        peptide_seqs,
        max_subs: int = 3,
        max_indels: int = 0,
        budget=None,
        threads: int = 0,
    ) -> pl.DataFrame:
        """Every alignment of every peptide within the substitution ceiling.

        Args:
            peptide_ids: Peptide ids, positionally aligned with ``peptide_seqs``.
            peptide_seqs: Peptide sequences, not yet collapsed.
            max_subs: Substitution ceiling. Used for all peptides unless
                ``budget`` is given.
            max_indels: Optional insertion/deletion ceiling. Gapped searches
                expose native reference spans and ``n_ins``/``n_dels`` counts;
                placements with indels are search-only until per-column
                alignment geometry is available for residue scoring.
            budget: Optional ``len -> max_subs`` callable, applied per peptide.
                Peptides are grouped by budget and each group is one batched
                call. Pass :func:`shipped_budget` to reproduce the upstream join.
            threads: 0 uses the process CPU budget; larger requests are clamped.

        Returns:
            One row per ``(clone_id, peptide_id, start, end, n_subs)`` alignment,
            plus ``mismatch_positions`` relative to the peptide. Peptides shorter than 6 residues are dropped -- the
            index cannot answer them -- and so are any the caller passes with
            characters outside the amino-acid alphabet.
        """
        threads = bounded_threads(threads)
        ids = list(peptide_ids)
        seqs = ms_collapse(peptide_seqs)
        if len(ids) != len(seqs):
            raise ValueError("peptide IDs and sequences must have the same length")
        keep = [i for i, s in enumerate(seqs) if len(s) >= MIN_QUERY_LEN]

        groups: dict[int, list[int]] = {}
        for i in keep:
            b = max_subs if budget is None else budget(len(seqs[i]))
            groups.setdefault(b, []).append(i)

        frames = []
        schema = {
            "clone_id": pl.String, "peptide_id": pl.String,
            "start": pl.Int32, "end": pl.Int32, "n_subs": pl.Int16,
            "mismatch_positions": pl.List(pl.Int16),
        }
        if max_indels:
            schema.update(n_ins=pl.Int32, n_dels=pl.Int32)
        clone_ids = pl.Series(self.clone_ids, dtype=pl.String)
        for b, idxs in sorted(groups.items()):
            res = self._index.search_batch(
                [seqs[i] for i in idxs], max_subs=b,
                max_indels=max_indels, threads=threads,
            )
            # seqtree 1.1 holds hits and mismatch coordinates in flat arrays.
            # Preserve its ordering without constructing a Python object per hit.
            arrays = res.to_numpy()
            query_row = np.repeat(np.arange(len(idxs)),
                                  np.diff(arrays["query_begin"]))
            # Widen before adding lengths; the strict output schema rejects
            # out-of-range coordinates/counts instead of NumPy wrapping them.
            start = arrays["offset"].astype(np.int64)
            lengths = arrays["length"].astype(np.int64)
            mismatches = pa.ListArray.from_arrays(
                pa.array(arrays["mm_begin"], type=pa.int32()),
                pa.array(arrays["mm_pos"], type=pa.int16()),
            )
            frames.append(pl.DataFrame({
                "clone_id": clone_ids.gather(arrays["ref_id"]),
                "peptide_id": pl.Series([ids[i] for i in idxs], dtype=pl.String).gather(query_row),
                "start": start,
                "end": start + lengths,
                "n_subs": arrays["n_subs"],
                "mismatch_positions": pl.from_arrow(mismatches),
                **({name: arrays[name] for name in ("n_ins", "n_dels")} if max_indels else {}),
            }, schema=schema))
        return pl.concat(frames) if frames else pl.DataFrame(schema=schema)
