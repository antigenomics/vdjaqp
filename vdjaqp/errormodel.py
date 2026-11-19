"""The de novo sequencing error model."""

from __future__ import annotations

import numpy as np
import polars as pl

from .alphabet import AA, N_AA, encode, encode_frame
from .data import load_clonotypes, load_peptides, load_targets, ms_stem_of
from .match import PeptideMatcher

__all__ = ["RESIDUE_MASS", "mass_prior", "isotonic", "ErrorModel", "fit_rna_junction_reliability"]

#: Monoisotopic residue masses, Da. Isoleucine is absent: it is merged into
#: leucine by the alphabet, being exactly isobaric with it.
RESIDUE_MASS = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919, "L": 113.08406, "N": 114.04293,
    "D": 115.02694, "Q": 128.05858, "K": 128.09496, "E": 129.04259,
    "M": 131.04049, "H": 137.05891, "F": 147.06841, "R": 156.10111,
    "Y": 163.06333, "W": 186.07931,
}

#: Mass scale of the decay, Da. Pairs closer than this are near-isobaric and
#: are the confusions de novo sequencing actually makes.
_MASS_SCALE = 0.5

#: Enrichment of near-isobaric confusions over what composition alone predicts.
#: Set so the prior reproduces the measured rate: 20.8% of observed
#: substitutions are under 1 Da, against 3.2% for a random amino-acid pair.
_MASS_ENRICHMENT = 12.0

#: Error rate for ``db`` rows, whose ``local_confidence`` is a constant 1.0 and
#: therefore carries no information. Fit by :meth:`ErrorModel.calibrate`.
_E_DB_DEFAULT = 0.10

#: Weight of the mass prior against empirical counts, in pseudo-observations.
_KAPPA = 50.0

#: Error rate is clamped here: no residue is certain, none is hopeless.
_E_MIN, _E_MAX = 0.005, 0.60


def mass_prior(composition: np.ndarray | None = None) -> np.ndarray:
    """``(19, 19)`` row-stochastic confusion prior, zero on the diagonal.

    ``P(called b | true a, an error occurred)``, as a *mixture*: a baseline that
    follows amino-acid composition, plus an excess concentrated on near-isobaric
    pairs. The mixture is what makes this usable as a prior -- a bare
    ``exp(-|dm| / 0.5)`` underflows to 1e-100 for distant pairs, so the cells the
    supervised substitutions never visit would cost hundreds of nats each and a
    single ordinary mismatch could veto a clonotype. De novo sequencing does
    make distant errors: the measured median ``|dm|`` is 16 Da.

    Args:
        composition: Amino-acid frequencies over :data:`AA`. Defaults to uniform.
    """
    if composition is None:
        composition = np.full(N_AA, 1.0 / N_AA)
    masses = np.array([RESIDUE_MASS[a] for a in AA])
    delta = np.abs(masses[:, None] - masses[None, :])
    matrix = composition[None, :] * (
        1.0 + _MASS_ENRICHMENT * np.exp(-delta / _MASS_SCALE)
    )
    np.fill_diagonal(matrix, 0.0)
    return matrix / matrix.sum(axis=1, keepdims=True)


def isotonic(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Monotone non-increasing fit of ``y`` on ``x`` by pool-adjacent-violators.

    Returns the sorted ``x`` breakpoints and the fitted values. Used to map
    de novo confidence to an error rate without assuming a functional form:
    error must fall as confidence rises, and nothing more is imposed.
    """
    order = np.argsort(x, kind="stable")
    xs, ys = x[order], y[order].astype(np.float64)
    values = list(ys)
    weights = [1.0] * len(ys)
    level = []
    for value, weight in zip(values, weights):
        level.append([value, weight])
        while len(level) > 1 and level[-2][0] < level[-1][0]:
            v2, w2 = level.pop()
            v1, w1 = level.pop()
            level.append([(v1 * w1 + v2 * w2) / (w1 + w2), w1 + w2])
    fitted = np.concatenate([np.full(int(w), v) for v, w in level])
    return xs, fitted


class ErrorModel:
    """``P(called residue | true residue, confidence)``.

    Build with :meth:`calibrate`. ``exclude`` holds out antibodies so a target
    is never scored with an error model its own peptides helped fit.

    Args:
        conf_x, conf_e: knots of the per-residue error rate as a function of de novo
            residue confidence (linear interpolation, constant beyond the end knots).
        substitution: substitution matrix, e.g. :func:`mass_prior`.
        n_substitutions: number of substitutions the calibration observed; recorded,
            not used in scoring.
        e_db: per-residue error rate for database-search hits, which carry no residue
            confidence. Pass it by keyword; a rate passed fourth positionally lands in
            ``n_substitutions`` and leaves ``e_db`` at its default.
    """

    def __init__(self, conf_x: np.ndarray, conf_e: np.ndarray, substitution: np.ndarray,
                 n_substitutions: int, e_db: float = _E_DB_DEFAULT):
        self._conf_x = conf_x
        self._conf_e = conf_e
        self.substitution = substitution
        self.n_substitutions = n_substitutions
        self.e_db = e_db

    # -- calibration -----------------------------------------------------
    @classmethod
    def calibrate(cls, exclude: set[str] | None = None, max_subs: int = 3,
                  root=None) -> "ErrorModel":
        """Fit on ``dataset_single_mabs`` self-matches.

        A peptide from antibody *X*'s own run, aligned to antibody *X*'s own
        clonotype, is a supervised observation of the error process: we know
        what the true residue was. ``db`` rows are excluded -- their
        disagreements measure the search software, not the instrument, and
        their ``alc`` is a constant 1.0.

        Args:
            exclude: ``clone_id`` values to leave out, for leave-one-antibody-out
                evaluation.
            max_subs: Search ceiling. Wider than the shipped ``floor(L/6)`` so
                the budget does not select away high-error peptides and bias
                the rate downward.
        """
        exclude = exclude or set()
        clonotypes = load_clonotypes("dataset_single_mabs", root)
        peptides = load_peptides("dataset_single_mabs", root)
        targets = load_targets(root)
        stem_to_clone = dict(zip(targets["ms_stem"], targets["clone_id"]))

        peptides = peptides.with_columns(ms_stem_of(peptides))
        peptides = peptides.with_columns(
            true_clone=pl.col("ms_stem").replace_strict(stem_to_clone, default=None)
        )
        known = peptides.filter(
            pl.col("true_clone").is_not_null()
            & ~pl.col("true_clone").is_in(list(exclude))
        )
        usable = known.filter(pl.col("search_type") == "denovo")

        # Only self-reference hits survive below. Keep the full clonotype table
        # for the composition prior, but do not search its discarded background.
        calibration_texts = clonotypes.filter(
            pl.col("clone_id").is_in(known["true_clone"].unique().implode())
        )
        matcher = PeptideMatcher.from_clonotypes(calibration_texts)

        def _self_hits(frame):
            got = matcher.search(frame["peptide_id"], frame["seq"], max_subs=max_subs)
            return got.join(
                frame.select(["peptide_id", "true_clone", "seq", "local_confidence"]),
                on="peptide_id",
            ).filter(pl.col("clone_id") == pl.col("true_clone"))

        hits = _self_hits(usable)
        # Keep the best offset per (peptide, clonotype): the true placement.
        # maintain_order on both, or the surviving row depends on thread
        # scheduling and the fitted error model stops being reproducible.
        hits = hits.sort("n_subs", maintain_order=True).unique(
            subset=["peptide_id", "clone_id"], keep="first", maintain_order=True
        )

        clone_code, clone_offset, _ = encode_frame(clonotypes["seq"])
        clone_row = (
            clonotypes.select("clone_id").with_row_index("clone_row")
        )
        hits = hits.join(clone_row, on="clone_id", how="left")

        pep_code, pep_offset, pep_length = encode_frame(hits["seq"])
        conf = np.concatenate(
            [np.asarray(x, dtype=np.float64) for x in hits["local_confidence"]]
        )
        starts = hits["start"].to_numpy().astype(np.int64)
        rows = hits["clone_row"].to_numpy().astype(np.int64)

        alignment_of = np.repeat(np.arange(hits.height), pep_length)
        ref_slot = (
            clone_offset[rows][alignment_of]
            + starts[alignment_of]
            + (np.arange(int(pep_length.sum()))
               - np.repeat(pep_offset, pep_length))
        )
        truth = clone_code[ref_slot]
        called = pep_code
        valid = (truth >= 0) & (called >= 0)
        wrong = valid & (truth != called)

        counts = np.zeros((N_AA, N_AA), dtype=np.float64)
        np.add.at(counts, (truth[wrong], called[wrong]), 1.0)
        conf = conf[valid]
        wrong = wrong[valid]

        conf_x, conf_e = isotonic(conf, wrong.astype(np.float64))
        conf_e = np.clip(conf_e, _E_MIN, _E_MAX)

        np.fill_diagonal(counts, 0.0)
        composition = np.bincount(
            clone_code[clone_code >= 0], minlength=N_AA
        ).astype(np.float64)
        composition /= composition.sum()

        prior = mass_prior(composition)
        substitution = (counts + _KAPPA * prior)
        substitution /= substitution.sum(axis=1, keepdims=True)

        # `db` rows report alc == 1.0 regardless, so their confidence carries no
        # information and a flat rate is fitted from their own self-matches.
        db = known.filter(pl.col("search_type") == "db")
        e_db = _E_DB_DEFAULT
        if db.height:
            db_hits = _self_hits(db)
            residues = sum(len(s) for s in db_hits["seq"])
            if residues:
                e_db = float(np.clip(
                    db_hits["n_subs"].sum() / residues, _E_MIN, _E_MAX
                ))
        return cls(conf_x, conf_e, substitution, int(counts.sum()), e_db)

    # -- query -----------------------------------------------------------
    def error_rate(self, confidence: np.ndarray, is_db: np.ndarray | None = None
                   ) -> np.ndarray:
        """Calibrated ``P(this residue was called wrong)``.

        Args:
            confidence: Per-residue ``local_confidence``.
            is_db: Mask of residues from ``search_type == 'db'`` peptides, whose
                confidence is a constant 1.0 and is replaced by the flat
                :attr:`e_db`. Leaving these at face value makes a database hit
                look infallible and turns each of its mismatches into a veto.
        """
        rate = np.interp(np.asarray(confidence, dtype=np.float64),
                         self._conf_x, self._conf_e)
        if is_db is not None:
            rate = np.where(is_db, self.e_db, rate)
        return rate

    def emission(self, confidence: np.ndarray) -> np.ndarray:
        """``(n, 19, 19)`` ``P(called | true)`` for each confidence value.

        Diagonal is ``1 - e``; the off-diagonal mass ``e`` is spread by the
        substitution matrix.
        """
        e = self.error_rate(confidence)[:, None, None]
        eye = np.eye(N_AA)[None, :, :]
        return (1.0 - e) * eye + e * self.substitution[None, :, :]


def fit_rna_junction_reliability(support, genuine) -> tuple[np.ndarray, np.ndarray]:
    """Optional monotone RNA-junction validity calibration, one row per class.

    Labels concern all generating RNA classes, never physical protein membership.
    For each tied nonnegative integer support, smooth the error fraction as
    (errors + .5)/(classes + 1). Repeat that mean once per class before the
    existing non-increasing PAV, preserving class-count weights and tie equality.
    Returned unique knots give genuine-class reliability strictly between 0/1.
    """
    x=np.asarray(support);y=np.asarray(genuine)
    if x.ndim!=1 or y.ndim!=1 or not x.size or x.size!=y.size:
        raise ValueError('support and genuine must be equal nonempty one-dimensional arrays')
    if x.dtype.kind not in 'iuf' or y.dtype.kind not in 'biuf':
        raise ValueError('support must be integer counts and genuine must be binary labels')
    if not np.isfinite(x).all() or (x<0).any() or (x>2**53).any() or (x!=np.floor(x)).any():
        raise ValueError('support must be finite nonnegative exactly representable integer counts')
    if not np.isfinite(y).all() or not np.isin(y,[0,1]).all():
        raise ValueError('genuine must contain binary RNA-class labels')
    knots,inverse,counts=np.unique(x.astype(np.float64),return_inverse=True,return_counts=True)
    errors=np.bincount(inverse,weights=1.-y.astype(np.float64),minlength=knots.size)
    smoothed=(errors+.5)/(counts+1.)
    _,fitted=isotonic(np.repeat(knots,counts),np.repeat(smoothed,counts))
    starts=np.r_[0,np.cumsum(counts)[:-1]]
    q=1.-fitted[starts]
    if not np.isfinite(q).all() or ((q<=0)|(q>=1)).any():
        raise ValueError('RNA reliability calibration did not produce finite probabilities in (0,1)')
    return knots,q
