"""Which peptides carry junction-specific signal?"""

from __future__ import annotations

import random

import numpy as np
import polars as pl
import seqtree

from .match import MIN_QUERY_LEN, ms_collapse
from .resources import bounded_threads

__all__ = ["shuffle_junctions", "peptide_evalues"]

#: ``region_ids`` value for the non-templated junction core.
_NDN = 6


def shuffle_junctions(clonotypes: pl.DataFrame, seed: int = 0) -> list[str]:
    """Each clonotype's sequence with its NDN core permuted in place.

    Framework and the V/J-templated flanks are untouched, so the decoy has the
    same residue composition, the same record lengths and the same germline
    redundancy as the repertoire -- it differs only where clonotypes differ from
    each other.

    Args:
        clonotypes: Needs ``seq`` and ``region_ids``.
        seed: Fixed so the decoy, and therefore every E-value, is reproducible.
    """
    rng = random.Random(seed)
    out = []
    for seq, region in zip(clonotypes["seq"].to_list(),
                           clonotypes["region_ids"].to_list()):
        chars = list(seq)
        core = np.flatnonzero(np.asarray(region) == _NDN)
        if core.size > 1:
            block = [chars[i] for i in core]
            rng.shuffle(block)
            for position, residue in zip(core, block):
                chars[position] = residue
        out.append("".join(chars))
    return out


def peptide_evalues(clonotypes: pl.DataFrame, peptides: pl.DataFrame,
                    max_subs: int = 2, seed: int = 0, k: int = 4,
                    threads: int | None = None) -> pl.DataFrame:
    """E-value per peptide against a junction-shuffled decoy repertoire.

    Args:
        clonotypes: The repertoire.
        peptides: Needs ``peptide_id`` and ``seq``.
        max_subs: Substitution radius, applied identically to both indexes.
        seed: Decoy shuffle seed.
        k: Seed width of the text indexes.
        threads: Defaults to this process's CPU budget.

    Returns:
        ``peptide_id``, ``n_target``, ``n_control``, ``E``, ``p_any``,
        ``p_enrichment``, ``rule_of_three``. Counts are distinct clonotypes.
        Peptides shorter than the index can answer are returned with null
        statistics rather than dropped.
    """
    threads = bounded_threads(threads)
    target = seqtree.TextIndex.build(
        ms_collapse(clonotypes["seq"].to_list()), "aa", k=k
    )
    decoy = seqtree.TextIndex.build(
        ms_collapse(shuffle_junctions(clonotypes, seed=seed)), "aa", k=k
    )

    sequences = ms_collapse(peptides["seq"].to_list())
    answerable = [
        i for i, s in enumerate(sequences) if len(s) >= max(MIN_QUERY_LEN, k)
    ]
    results = seqtree.text_evalues(
        target, decoy, [sequences[i] for i in answerable],
        max_subs=max_subs, threads=threads,
    )

    fields = ("n_target", "n_control", "E", "p_any", "p_enrichment",
              "rule_of_three")
    columns = {name: [None] * peptides.height for name in fields}
    for position, record in zip(answerable, results):
        for name in fields:
            columns[name][position] = record[name]
    return peptides.select("peptide_id").with_columns(
        **{name: pl.Series(values) for name, values in columns.items()}
    )
