"""The scoring pipeline as one object: peptides in, scored clonotypes out."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import polars as pl

from .aggregate import ClonotypeScorer
from .background import Background
from .data import load_clonotypes, load_peptides
from .errormodel import ErrorModel
from .evidence import AlignmentScorer
from .match import PeptideMatcher, shipped_budget

__all__ = ["Evidence"]


@dataclass(frozen=True)
class Evidence:
    """Every intermediate of one dataset's scoring run.

    Attributes
    ----------
    dataset
        Which of :data:`vdjaqp.data.DATASETS` this was built from.
    clonotypes, peptides
        The two inputs, as loaded. ``peptides`` may be a subset of the
        dataset's own table -- the leave-one-out evaluation passes a single
        run's peptides here.
    matches
        One row per (peptide, clonotype, offset) alignment from the substring
        search, before any scoring.
    alignments
        ``matches`` with a per-alignment log likelihood ratio attached.
    background, error
        The fitted components, exposed because they are results in their own
        right: ``background`` carries the per-region information content and
        ``error`` the de novo substitution matrix.
    """

    dataset: str
    clonotypes: pl.DataFrame
    peptides: pl.DataFrame
    matches: pl.DataFrame
    alignments: pl.DataFrame
    background: Background
    error: ErrorModel

    @classmethod
    def build(cls, dataset: str = "dataset_mixa_gel", *, root=None,
              peptides: pl.DataFrame | None = None,
              max_subs: int | None = None,
              exclude: Iterable[str] = ()) -> "Evidence":
        """Run the search and the per-alignment scoring.

        Parameters
        ----------
        peptides
            Use this peptide table instead of the dataset's own. The
            leave-one-out evaluation passes one run's peptides.
        max_subs
            Fixed substitution budget. The default follows the upstream
            ``floor(L/6)`` rule.
        exclude
            Clone ids kept out of the de novo error-model calibration, so an
            antibody's own peptides never inform the model that scores it.
        """
        clonotypes = load_clonotypes(dataset, root)
        if peptides is None:
            peptides = load_peptides(dataset, root)
        error = ErrorModel.calibrate(exclude=set(exclude), root=root)
        return cls.from_frames(clonotypes, peptides, error_model=error,
                               max_subs=max_subs, dataset=dataset)

    @classmethod
    def from_frames(cls, clonotypes: pl.DataFrame, peptides: pl.DataFrame, *,
                    error_model: ErrorModel, max_subs: int | None = None,
                    dataset: str = "frames", include_components: bool = False) -> "Evidence":
        """Score supplied observations with an explicit, independently defined model.

        No loading/calibration or simulator truth inference occurs here. This
        shares the exact search/background/likelihood route with ``build``.
        ``include_components=True`` retains optional observed-AA emission and
        original context components in the same single alignment-scoring call.
        """
        matcher = PeptideMatcher.from_clonotypes(clonotypes)
        if max_subs is None:
            matches = matcher.search(peptides["peptide_id"], peptides["seq"],
                                     budget=shipped_budget)
        else:
            matches = matcher.search(peptides["peptide_id"], peptides["seq"],
                                     max_subs=max_subs)
        error = error_model
        background = Background(clonotypes)
        alignments = AlignmentScorer(clonotypes, error, background).score(
            matches, peptides, include_components=include_components)
        return cls(dataset=dataset, clonotypes=clonotypes, peptides=peptides,
                   matches=matches, alignments=alignments,
                   background=background, error=error)

    def score(self, **kw) -> pl.DataFrame:
        """Aggregate alignments to one score per clonotype.

        Keyword arguments go to :meth:`vdjaqp.aggregate.ClonotypeScorer.score`
        (``use_prior``, ``competition``), which is how the ablation separates
        the contribution of each component.
        """
        return ClonotypeScorer(self.clonotypes, self.background).score(
            self.alignments, self.peptides, **kw)

    def match_counts(self) -> pl.DataFrame:
        """The weak baseline: how many peptides matched, and nothing else.

        Clonotypes with no match are carried at 0 rather than dropped, so the
        baseline ranks the same candidate rows the model does.
        """
        counts = self.matches.group_by("clone_id").agg(
            score=pl.len().cast(pl.Float64))
        missing = self.clonotypes.select("clone_id").join(
            counts, on="clone_id", how="anti")
        return pl.concat([counts, missing.with_columns(score=pl.lit(0.0))])
