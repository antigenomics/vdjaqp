"""Expanding alignments to residue slots, without a Python loop per alignment."""

from __future__ import annotations

import numpy as np

__all__ = ["expand_alignments"]


def expand_alignments(pep_offset, pep_length, pep_row, clone_offset, clone_row,
                      start):
    """Residue-level index arrays for a block of alignments.

    Args:
        pep_offset: Start of each peptide's block in the flat peptide arrays.
        pep_length: Length of each peptide.
        pep_row: Peptide index of each alignment.
        clone_offset: Start of each clonotype's block in the flat clonotype
            arrays.
        clone_row: Clonotype index of each alignment.
        start: Offset within the clonotype at which each alignment begins.

    Returns:
        ``(alignment_of, pep_slot, ref_slot)``, each of length
        ``sum(pep_length[pep_row])``. ``alignment_of`` says which alignment a
        residue slot belongs to, so a reduction back to per-alignment values is
        one ``np.add.at``.
    """
    lengths = pep_length[pep_row]
    total = int(lengths.sum())
    alignment_of = np.repeat(np.arange(len(pep_row)), lengths)
    block_start = np.concatenate([[0], np.cumsum(lengths)[:-1]])
    within = np.arange(total) - np.repeat(block_start, lengths)
    pep_slot = pep_offset[pep_row][alignment_of] + within
    ref_slot = clone_offset[clone_row][alignment_of] + start[alignment_of] + within
    return alignment_of, pep_slot, ref_slot
