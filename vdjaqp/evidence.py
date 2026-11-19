"""Per-alignment log likelihood ratio."""

from __future__ import annotations

import numpy as np
import polars as pl

from .alphabet import encode_frame
from .background import Background
from .errormodel import ErrorModel
from .layout import expand_alignments

__all__ = ["AlignmentScorer"]

#: Alignments scored per chunk. Each alignment expands to one row per residue,
#: so a chunk of 2e6 alignments is about 3e7 residue slots and a few hundred
#: megabytes of scratch -- bounded regardless of how wide the search was.
_CHUNK = 2_000_000


class AlignmentScorer:
    """Scores every alignment in a match frame with its log likelihood ratio."""

    def __init__(self, clonotypes: pl.DataFrame, error_model: ErrorModel,
                 background: Background | None = None):
        self.error = error_model
        self.background = background or Background(clonotypes)
        self._clone_index = clonotypes.select("clone_id", pl.col("seq").str.len_chars().alias("_reference_length")).with_row_index(
            "clone_row"
        ).with_columns(pl.col("clone_row").cast(pl.Int64))

    def score(self, matches: pl.DataFrame, peptides: pl.DataFrame, *,
              include_components: bool = False) -> pl.DataFrame:
        """Add ``llr`` (nats) to a match frame.

        Args:
            matches: Output of :meth:`vdjaqp.match.PeptideMatcher.search`.
            peptides: Must carry ``peptide_id``, ``seq``, ``local_confidence``
                and ``search_type``.
            include_components: Also expose observed-AA ``log_emission`` and
                the original context-masked ``log_context_emission`` and
                ``log_context_null``. Noncanonical observed/source residues
                make the raw emission unavailable (NaN), rather than favoring
                an incomplete likelihood. Defaults and ``llr`` stay unchanged.

        Returns:
            ``matches`` with ``llr`` appended, row order preserved. Positive
            means the clonotype explains the reading better than the repertoire
            at large does.
        """
        if any(name in matches.columns and (matches[name].is_null() | matches[name].ne(0)).any()
               for name in ('n_ins', 'n_dels')):
            raise ValueError('gapped placements are search-only; residue scoring requires per-column alignment geometry')
        pep_code, pep_offset, pep_length = encode_frame(peptides["seq"])
        pep_conf = np.concatenate(
            [np.asarray(x, dtype=np.float32) for x in peptides["local_confidence"]]
        ) if peptides.height else np.empty(0, dtype=np.float32)
        pep_index = peptides.select("peptide_id", pl.col("seq").str.len_chars().alias("_peptide_length")).with_row_index("pep_row")
        is_db = (peptides["search_type"] == "db").to_numpy()

        keyed = (
            matches.with_row_index("_align_row")
            .join(self._clone_index, on="clone_id", how="left")
            .join(pep_index, on="peptide_id", how="left")
        )
        if keyed["clone_row"].null_count() or keyed["pep_row"].null_count():
            raise ValueError("match frame references unknown clone_id or peptide_id")
        start, end = pl.col('start').cast(pl.Int64), pl.col('end').cast(pl.Int64)
        if keyed.select(((start < 0) | (end < start) |
                         (end-start != pl.col('_peptide_length')) |
                         (end > pl.col('_reference_length'))).fill_null(True).any()).item():
            raise ValueError('invalid ungapped placement span or reference bounds')

        llr = np.empty(keyed.height, dtype=np.float64)
        components={name:np.empty(keyed.height,dtype=np.float64) for name in
            ('log_emission','log_context_emission','log_context_null')} if include_components else {}
        for begin in range(0, keyed.height, _CHUNK):
            chunk = keyed.slice(begin, _CHUNK)
            values = self._score_chunk(
                chunk, pep_code, pep_offset, pep_length, pep_conf, is_db,
                include_components=include_components)
            if include_components:
                values,extra=values
                for name,array in extra.items():components[name][begin:begin+chunk.height]=array
            llr[begin:begin + chunk.height] = values
        return keyed.sort("_align_row").drop(
            ["_align_row", "clone_row", "pep_row", "_reference_length", "_peptide_length"]
        ).with_columns(llr=pl.Series(llr, dtype=pl.Float64),
                       **{name:pl.Series(array,dtype=pl.Float64) for name,array in components.items()})

    def _score_chunk(self, chunk, pep_code, pep_offset, pep_length, pep_conf,
                     is_db, *, include_components=False):
        """Log likelihood ratio for one bounded block of alignments."""
        pep_row = chunk["pep_row"].to_numpy().astype(np.int64)
        clone_row = chunk["clone_row"].to_numpy()
        start = chunk["start"].to_numpy().astype(np.int64)

        alignment_of, pep_slot, ref_slot = expand_alignments(
            pep_offset, pep_length, pep_row,
            self.background.offset, clone_row, start,
        )

        called = pep_code[pep_slot].astype(np.int64)
        truth = self.background.code[ref_slot].astype(np.int64)
        valid = (called >= 0) & (truth >= 0) & self.background.context_valid[ref_slot]
        called_safe = np.where(valid, called, 0)
        truth_safe = np.where(valid, truth, 0)

        e = self.error.error_rate(
            pep_conf[pep_slot], is_db[pep_row][alignment_of]
        )
        numerator = (1.0 - e) * (called_safe == truth_safe) + e * (
            self.error.substitution[truth_safe, called_safe]
        )
        denominator = (
            (1.0 - e) * self.background.frequency(ref_slot, called_safe)
            + e * self.background.mixed(
                ref_slot, called_safe, self.error.substitution
            )
        )
        per_residue = np.where(
            valid, np.log(numerator) - np.log(denominator), 0.0
        )
        out = np.zeros(chunk.height)
        np.add.at(out, alignment_of, per_residue)
        if not include_components:return out
        # Source compatibility uses observed amino acids even when their gene
        # context is unavailable. Context validity belongs only to the original
        # information LLR and must not make a known mismatch uninformative.
        observed_valid=(called>=0)&(truth>=0)
        observed_called=np.where(observed_valid,called,0)
        observed_truth=np.where(observed_valid,truth,0)
        emission=(1.-e)*(observed_called==observed_truth)+e*self.error.substitution[observed_truth,observed_called]
        terms={
            'log_emission':np.where(observed_valid,np.log(emission),0.),
            'log_context_emission':np.where(valid,np.log(numerator),0.),
            'log_context_null':np.where(valid,np.log(denominator),0.)}
        extra={name:np.bincount(alignment_of,weights=value,minlength=chunk.height)
               for name,value in terms.items()}
        observed_count=np.bincount(alignment_of,weights=observed_valid,minlength=chunk.height)
        span_count=np.bincount(alignment_of,minlength=chunk.height)
        extra['log_emission'][observed_count!=span_count]=np.nan
        return out,extra
