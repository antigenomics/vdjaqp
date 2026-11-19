"""Reverting somatic hypermutation before scoring a junction."""

from __future__ import annotations

import numpy as np
import polars as pl

from .data import to_gene_level
from .germline import aligned_columns, canonical_j, canonical_v, common_reference

__all__ = ["revert_junctions"]

#: ``region_ids`` value for the non-templated junction core.
_NDN = 6

#: How far past each CDR3 boundary to align, so the templated flanks of the
#: junction are covered even when V or J contributes unusually many residues.
_V_WINDOW, _J_WINDOW = 8, 10


def revert_junctions(clonotypes: pl.DataFrame) -> pl.DataFrame:
    """Replace hypermutated residues in each junction's templated flanks.

    Args:
        clonotypes: Needs ``clone_id``, ``seq``, ``region_ids``, ``cdr3_start``,
            ``cdr3_end``, ``v_name`` and ``j_name``.

    Returns:
        ``clone_id``, ``junction`` (as observed), ``junction_reverted``,
        ``n_reverted`` (residues changed) and ``reverted`` (whether any were).
        Row order matches the input.
    """
    frame = to_gene_level(clonotypes)
    v_germline, j_germline = canonical_v(), canonical_j()

    observed, restored, changed = [], [], []
    for seq, region, start, end, v_gene, j_gene in zip(
        frame["seq"].to_list(), frame["region_ids"].to_list(),
        frame["cdr3_start"].to_numpy(), frame["cdr3_end"].to_numpy(),
        frame["v_name"].to_list(), frame["j_name"].to_list(),
    ):
        start, end = int(start), int(end)
        junction = seq[start:end]
        observed.append(junction)
        core = np.asarray(region)[start:end] == _NDN
        residues = list(junction)

        v_match = common_reference(v_gene, v_germline)
        if v_match is not None:
            _, v_reference = v_match
            stop = min(len(seq), start + _V_WINDOW)
            for position, _, germline_residue in aligned_columns(seq[:stop],
                                                                 v_reference):
                index = position - start
                if 0 <= index < len(residues) and not core[index]:
                    residues[index] = germline_residue

        j_match = common_reference(j_gene, j_germline)
        if j_match is not None:
            _, j_reference = j_match
            begin = max(0, end - _J_WINDOW)
            for position, _, germline_residue in aligned_columns(seq[begin:],
                                                                 j_reference):
                index = begin + position - start
                if 0 <= index < len(residues) and not core[index]:
                    residues[index] = germline_residue

        reverted = "".join(residues)
        restored.append(reverted)
        changed.append(sum(a != b for a, b in zip(junction, reverted)))

    return frame.select("clone_id").with_columns(
        junction=pl.Series(observed),
        junction_reverted=pl.Series(restored),
        n_reverted=pl.Series(changed, dtype=pl.Int32),
        reverted=pl.Series(np.asarray(changed) > 0),
    )
