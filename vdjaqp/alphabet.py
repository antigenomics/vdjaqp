"""The mass-spectrometry amino-acid alphabet, encoded in whole batches."""

from __future__ import annotations

import numpy as np
import polars as pl

__all__ = ["AA", "N_AA", "encode", "encode_frame"]

#: The 19 distinct residue masses. Isoleucine is folded into leucine.
AA = "ACDEFGHKLMNPQRSTVWY"
N_AA = len(AA)

#: Byte value -> alphabet index, with -1 for anything outside it (``X``, ``*``,
#: ``B``, ``Z`` and lower case all land here and are skipped, not guessed).
_LUT = np.full(256, -1, dtype=np.int8)
for _i, _a in enumerate(AA):
    _LUT[ord(_a)] = _i
_LUT[ord("I")] = _LUT[ord("L")]


def encode(seq: str) -> np.ndarray:
    """Encode one sequence to alphabet indices."""
    return _LUT[np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)]


def encode_frame(seqs: pl.Series) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Encode a whole column of sequences in one pass.

    Args:
        seqs: Amino-acid strings.

    Returns:
        ``(code, offset, length)``. ``code`` is the concatenation of every
        sequence's indices; sequence ``i`` occupies
        ``code[offset[i]:offset[i] + length[i]]``.
    """
    length = seqs.str.len_bytes().to_numpy().astype(np.int64)
    offset = np.concatenate([[0], np.cumsum(length)[:-1]]) if len(length) else length
    joined = "".join(seqs.to_list())
    code = _LUT[np.frombuffer(joined.encode("ascii", "replace"), dtype=np.uint8)]
    return code, offset, length
