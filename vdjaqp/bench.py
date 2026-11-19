"""Where do the 12 known antibodies land, and what did each component buy?"""

from __future__ import annotations

import numpy as np
import polars as pl

from .data import load_peptides, load_targets, ms_stem_of, to_gene_level
from .pipeline import Evidence

__all__ = ["classification_metrics", "rank_targets", "summarise", "ablation", "leave_one_mab_out"]


def rank_targets(table: pl.DataFrame, score_column: str, targets: set[str]
                 ) -> pl.DataFrame:
    """Rank of every target in a scored clonotype table.

    Ties break on ``clone_id``, including for the baselines, so every method in
    a comparison is ranked the same way and a reported rank does not move
    between identical runs. ``n_tied`` says how many clonotypes share the
    score, because a rank drawn from inside a large tie block is arbitrary.
    """
    ranked = (
        table.sort([score_column, "clone_id"], descending=[True, False])
        .with_columns(n_tied=pl.len().over(score_column))
        .with_row_index("rank", offset=1)
    )
    return ranked.filter(pl.col("clone_id").is_in(targets)).sort("rank")


def summarise(ranks: np.ndarray, n: int) -> dict:
    """Median, worst and recall-at-k for one method."""
    return {
        "n_clonotypes": n,
        "median_rank": float(np.median(ranks)),
        "worst_rank": int(ranks.max()),
        "top_50": int((ranks <= 50).sum()),
        "top_500": int((ranks <= 500).sum()),
    }


def classification_metrics(table: pl.DataFrame, targets: set[str], score_column: str = "score") -> dict:
    """Record-label benchmark metrics; AP uses tied score groups, not trapezoidal PR area.

    Other records are benchmark negatives, not established biological nonbinders.
    Binary calls use the existing confirmed_nats >= 20 rule, independently of rank score.
    """
    if (not targets or not targets <= set(table["clone_id"])
            or table["clone_id"].null_count() or table["clone_id"].n_unique()!=table.height):
        raise ValueError("targets must be present and candidate IDs unique/nonmissing")
    score=table[score_column].to_numpy().astype(float)
    confirmed=table["confirmed_nats"].to_numpy().astype(float)
    if not np.isfinite(score).all() or not np.isfinite(confirmed).all():
        raise ValueError("benchmark scores must be finite")
    positive=table["clone_id"].is_in(targets).to_numpy()
    order=np.argsort(-score,kind="stable");score=score[order];y=positive[order].astype(int)
    ends=np.r_[np.flatnonzero(score[1:]!=score[:-1])+1,len(score)]
    tp=np.cumsum(y)[ends-1];fp=ends-tp;p=len(targets);n=table.height-p
    recall=tp/p;precision=tp/ends
    ap=float(np.sum(np.diff(np.r_[0.,recall])*precision))
    auc=float(np.sum(np.diff(np.r_[0.,fp/n])*(np.r_[0.,recall][:-1]+recall)/2)) if n else None
    calls=confirmed>=20;true=int(np.sum(calls & positive));false=int(np.sum(calls & ~positive));miss=p-true
    ranks=rank_targets(table,score_column,targets)["rank"].to_numpy()
    return dict(n_candidates=table.height,n_positives=p,n_benchmark_negatives=n,
        average_precision=ap,roc_auc=auc,threshold_confirmed_nats=20,
        true_positive=true,false_positive=false,false_negative=miss,true_negative=n-false,
        precision_at20=true/(true+false) if true+false else 0.,recall_at20=true/p,
        f1_at20=2*true/(2*true+false+miss),**summarise(ranks,table.height))


def ablation(dataset: str = "dataset_mixa_gel", max_subs: int | None = None,
             root=None) -> pl.DataFrame:
    """One row per method increment, so each component's effect is separable."""
    ev = Evidence.build(dataset, root=root, max_subs=max_subs)
    targets = set(load_targets(root)["clone_id"])

    rows = []
    counts = ev.match_counts()
    rows.append(("match count (baseline)", counts, "score"))

    isotype = to_gene_level(ev.clonotypes).select(
        ["clone_id", pl.col("c_name").str.starts_with("IGHG").alias("igg")]
    )
    rows.append((
        "match count + isotype",
        counts.join(isotype, on="clone_id").with_columns(
            score=pl.col("score") + pl.col("igg").cast(pl.Float64) * 1e6
        ),
        "score",
    ))
    rows.append((
        "LLR, no competition",
        ev.score(use_prior=False, competition=False), "score",
    ))
    rows.append((
        "LLR + competition",
        ev.score(use_prior=False), "score",
    ))
    rows.append((
        "LLR + competition + isotype",
        ev.score(use_prior=True), "score",
    ))

    out = []
    for name, table, column in rows:
        ranks = rank_targets(table, column, targets)["rank"].to_numpy()
        out.append({"method": name, **summarise(ranks, table.height),
                    "ranks": " ".join(str(int(r)) for r in ranks)})
    return pl.DataFrame(out)


def leave_one_mab_out(root=None) -> pl.DataFrame:
    """Per-antibody evaluation on the individual monoclonal runs.

    For each of the eight antibodies with its own MS run, the clonotypes are
    ranked using **only that run's peptides**, with the de novo error model
    refitted on the other seven. No antibody's own peptides ever inform the
    error model that scores it.

    ``dataset_single_mabs`` must be evaluated this way. Pooled it is
    meaningless -- every run's peptides are in the same table -- and
    ``source_file`` names the antibody, so it is used here to split runs and
    is never visible to a score.
    """
    peptides = load_peptides("dataset_single_mabs", root)
    single = load_targets(root).filter(pl.col("has_single_run") == "yes")
    peptides = peptides.with_columns(ms_stem_of(peptides))

    rows = []
    for stem, clone_id in zip(single["ms_stem"], single["clone_id"]):
        own = peptides.filter(pl.col("ms_stem") == stem).drop("ms_stem")
        ev = Evidence.build("dataset_single_mabs", root=root, peptides=own,
                            exclude={clone_id})
        table = ev.score()
        # Tie-break on clone_id, as rank_targets() does. Without it a rank drawn from
        # inside a tie block is arbitrary and moves between identical runs -- the
        # match-count baseline ties thousands of clonotypes at the same integer.
        ranked = table.sort(["score", "clone_id"], descending=[True, False]).with_row_index(
            "rank", offset=1)
        hit = ranked.filter(pl.col("clone_id") == clone_id)
        base = ev.match_counts().sort(
            ["score", "clone_id"], descending=[True, False]
        ).with_row_index("rank", offset=1)
        rows.append({
            "ms_stem": stem,
            "clone_id": clone_id,
            "n_peptides": own.height,
            "rank": int(hit["rank"][0]),
            "score_nats": float(hit["score"][0]),
            "rank_match_count": int(
                base.filter(pl.col("clone_id") == clone_id)["rank"][0]
            ),
        })
    return pl.DataFrame(rows).sort("rank")
