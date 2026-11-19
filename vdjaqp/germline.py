"""Canonical germline V and J amino-acid sequences, one per gene."""

from __future__ import annotations

from functools import lru_cache

import seqtree
from arda.germline import segment_germlines
from arda.refbuild.translate import translate
from seqtree.pairwise import align as _align

__all__ = ["canonical_v", "canonical_j", "GAP_OPEN", "GAP_EXTEND",
           "blosum62", "aligned_columns"]

#: Affine gap costs for clonotype-to-germline alignment, on BLOSUM62.
GAP_OPEN, GAP_EXTEND = 11, 1

#: Motif marking the J frame: the FR4 tryptophan and the GxG that follows it.
_J_MOTIF = "WG"


def common_reference(call: str | None, references: dict[str, str]) -> tuple[str, str] | None:
    """Resolve alternatives only when every gene has the same AA reference.

    The returned key is deterministic; the original producer call is untouched.
    Missing or sequence-distinct alternatives never select an arbitrary template.
    """
    genes = sorted({part.strip().split("*")[0] for part in (call or "").split(",")})
    sequences = [references.get(gene) for gene in genes]
    if sequences[0] is None or any(seq != sequences[0] for seq in sequences):
        return None
    return ",".join(genes), sequences[0]


def _gene_of(allele: str) -> str:
    return allele.split("*")[0]


def _pick_canonical(alleles: dict[str, str], locus: str) -> dict[str, list[str]]:
    """Group alleles by gene, lowest allele number first."""
    by_gene: dict[str, list[tuple[str, str]]] = {}
    for allele, seq in alleles.items():
        if allele.startswith(locus):
            by_gene.setdefault(_gene_of(allele), []).append((allele, seq))
    return {g: [s for _, s in sorted(v)] for g, v in by_gene.items()}


@lru_cache(maxsize=4)
def canonical_v(locus: str = "IGHV") -> dict[str, str]:
    """``{gene: germline amino-acid sequence}`` for the V segment.

    The scaffold reference is trimmed to the coding frame, so frame 0 is
    correct and the sequence runs from FR1 through the conserved Cys104.
    """
    out = {}
    for gene, seqs in _pick_canonical(segment_germlines("human", "v"), locus).items():
        out[gene] = translate(seqs[0])
    return out


@lru_cache(maxsize=4)
def canonical_j(locus: str = "IGHJ") -> dict[str, str]:
    """``{gene: germline amino-acid sequence}`` for the J segment.

    J scaffolds are not frame-trimmed, so the reading frame is chosen as the one
    carrying the FR4 ``WG`` motif -- translating frame 0 blindly yields nonsense
    like ``TTGSTPGAREPWSPSPQ`` for ``IGHJ5*02``.
    """
    out = {}
    for gene, seqs in _pick_canonical(segment_germlines("human", "j"), locus).items():
        nt = seqs[0]
        best = max(
            (translate(nt[f:]) for f in range(3)),
            key=lambda aa: (_J_MOTIF in aa, -aa.count("*")),
        )
        out[gene] = best
    return out


@lru_cache(maxsize=1)
def blosum62():
    """The substitution matrix every clonotype-to-germline alignment uses."""
    return seqtree.SubstitutionMatrix.blosum62()


def aligned_columns(query: str, reference: str):
    """Columns where a query residue faces a germline residue, no gap on either side.

    Yields ``(query_position, reference_position, reference_residue)``. This is
    the one place a clonotype is put into register with its germline: the
    background keys its frequency cells on ``reference_position`` and the SHM
    reversion substitutes ``reference_residue``, so both see exactly the same
    alignment and cannot drift apart.
    """
    if not query or not reference:
        return
    alignment = _align(query, reference, blosum62(), mode="global",
                       gap_open=GAP_OPEN, gap_extend=GAP_EXTEND)
    query_at = reference_at = 0
    for query_char, reference_char in zip(alignment.aligned_query,
                                          alignment.aligned_ref):
        if query_char != "-" and reference_char != "-":
            yield query_at, reference_at, reference_char
        query_at += query_char != "-"
        reference_at += reference_char != "-"
