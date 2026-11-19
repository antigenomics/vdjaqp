"""From per-alignment evidence to a ranked list of clonotypes."""

from __future__ import annotations

import numpy as np
import polars as pl

from .alphabet import encode_frame
from .background import Background
from .data import isotype_of, to_gene_level
from .layout import expand_alignments

__all__ = ["spectrum_key", "ClonotypeScorer", "compete_identity_classes",
           "fit_identity_mixture", "profile_identity_removals",
           "compete_identity_emissions", "identity_rna_support",
           "apply_rna_junction_reliability"]

#: Retention-time window, in minutes, for splitting a (file, mass, charge)
#: group into separate precursors. Matches the upstream build's own ``rt_tol``.
_RT_TOL = 2.0

#: Prior probability that a purified serum antibody is IgG. The repertoire is
#: 87% IgM, so this is the most informative non-sequence evidence available.
_P_IGG = 0.90

#: Caps, in nats. A confirmed residue is worth at most ``_W_MAX``; a
#: contradicted clonotype loses at most ``_D_MAX``. Both exist so that one long
#: framework peptide cannot decide a ranking of 33,052 clonotypes.
_W_MAX, _D_MAX = 8.0, 16.0

#: Confidence assumed for a peptide whose ``alc`` is missing. Database-search
#: rows carry no usable confidence -- the README describes them as a constant
#: 1.0, but they are NaN in two of the three datasets -- so they are given a
#: neutral value rather than being believed or dropped.
_ALC_DEFAULT = 0.5

#: Alignments per chunk, bounding the residue-level scratch arrays.
_CHUNK = 2_000_000

#: ``region_ids`` values that count as junctional: the templated CDR3 flanks and
#: the non-templated core. Confirming a residue here is evidence about *this*
#: clonotype; confirming framework is mostly evidence about its V gene.
_JUNCTION_REGIONS = (5, 6)

#: Floor, in nats, below which ``framework_only`` stays false. One nat is a
#: factor of e in likelihood; under that there is no claim to qualify, and
#: without a floor the flag fires on float dust -- a clonotype credited
#: 4e-11 nats has no junctional evidence only because it has none at all.
_FLAG_MIN_NATS = 1.0


def isotype_log_prior(isotypes) -> np.ndarray:
    """Retain the legacy complete-data prior; missing calls supply no evidence.

    With missing calls, known categories receive log(target/source frequency),
    estimated over observed calls only; unknown calls receive zero. This assumes
    observed isotypes describe the available category mix, not that missing calls
    are non-IgG. A single observed category cannot define relative odds.
    """
    isotypes = np.asarray(isotypes)
    known = isotypes != "unknown"
    is_igg = isotypes == "IgG"
    n_known, n_igg = int(known.sum()), int(is_igg.sum())
    if n_igg == 0 or n_igg == n_known:
        return np.zeros(len(isotypes))
    if known.all():
        log_odds = np.log((_P_IGG / n_igg) / ((1 - _P_IGG) / (len(is_igg) - n_igg)))
        return np.where(is_igg, log_odds, 0.0)
    prior = np.zeros(len(isotypes))
    prior[is_igg] = np.log(_P_IGG / (n_igg / n_known))
    prior[known & ~is_igg] = np.log((1 - _P_IGG) / ((n_known - n_igg) / n_known))
    return prior


def _rank_order(table: pl.DataFrame) -> pl.DataFrame:
    """Sort by score, breaking ties on ``clone_id`` so ranks are reproducible.

    Most of the repertoire ends on exactly the same score -- the isotype prior
    with no confirmed evidence -- and polars does not promise an order within
    ties. Without a tiebreaker a target's reported rank moved by several places
    between identical runs. ``n_tied`` is carried alongside so a rank inside a
    large tie block is not read as a measurement.
    """
    return (
        table.sort(["score", "clone_id"], descending=[True, False])
        .with_columns(n_tied=pl.len().over("score"))
    )


def spectrum_key(peptides: pl.DataFrame) -> pl.DataFrame:
    """Group readings that came from the same physical measurement.

    ``peptide_id`` is seven ``/``-joined fields: source file, mass, charge,
    m/z, retention time, scan, sequence. Readings of one precursor share file,
    mass and charge; retention time separates species that merely collide in
    mass and charge.
    """
    if "spectrum_id" in peptides.columns:
        if peptides.schema["spectrum_id"] != pl.String or peptides.select(
            (pl.col("spectrum_id").str.strip_chars() == "").fill_null(True).any()
        ).item():
            raise ValueError("spectrum_id must contain nonempty native observation IDs")
        return peptides
    fields = pl.col("peptide_id").str.split("/")
    return peptides.with_columns(
        spectrum_id=pl.concat_str(
            [
                fields.list.get(0),
                fields.list.get(1),
                fields.list.get(2),
                (fields.list.get(4).cast(pl.Float64, strict=False) / _RT_TOL)
                .floor()
                .cast(pl.Int64)
                .cast(pl.String),
            ],
            separator="|",
            ignore_nulls=True,
        )
    )


class ClonotypeScorer:
    """Turns scored alignments into a ranked clonotype table."""

    def __init__(self, clonotypes: pl.DataFrame, background: Background):
        self.clonotypes = to_gene_level(clonotypes)
        self.background = background
        self._clone_index = self.clonotypes.select("clone_id").with_row_index(
            "clone_row"
        ).with_columns(pl.col("clone_row").cast(pl.Int64))
        self._information = np.minimum(background.information(), _W_MAX)
        region = np.concatenate(
            [np.asarray(r) for r in self.clonotypes["region_ids"].to_list()]
        )
        self._is_junction = np.isin(region, _JUNCTION_REGIONS)

    def score(self, scored: pl.DataFrame, peptides: pl.DataFrame,
              use_prior: bool = True, competition: bool = True,
              aggregation: str = "alignment_union") -> pl.DataFrame:
        """Rank clonotypes by confirmed information content.

        Args:
            scored: Alignments carrying ``llr``, from
                :meth:`vdjaqp.evidence.AlignmentScorer.score`.
            peptides: The peptide table, for the spectrum grouping.
            use_prior: Add the isotype log-prior.
            competition: Resolve the spectrum-level competition. Turning it off
                is the ablation that shows what the competition is worth.
            aggregation: ``alignment_union`` preserves the historical heuristic;
                ``exclusive_spectrum`` sums agreeing alternative responsibilities
                within a spectrum before unioning across spectra. With competition
                disabled, both use the historical deterministic coverage union.

        Returns:
            One row per clonotype, sorted by ``score`` descending:

            * ``confirmed_nats`` -- information the mass spec confirmed.
            * ``junction_nats`` -- the part of it inside CDR3 or the
              non-templated core. This is the part that is about *this*
              clonotype rather than about its V gene.
            * ``contradiction``, ``evidence``, ``score`` -- in nats.
            * ``framework_only`` -- ``confirmed_nats`` at or above
              :data:`_FLAG_MIN_NATS` with ``junction_nats`` exactly zero. A true
              value means every nat this clonotype was credited came from
              residues it shares with the rest of its V gene, so the score must
              not be read as identification. The row is flagged, not removed.
        """
        if aggregation not in ("alignment_union", "exclusive_spectrum"):
            raise ValueError("aggregation must be alignment_union or exclusive_spectrum")
        peptides = spectrum_key(peptides)
        pep_code, pep_offset, pep_length = encode_frame(peptides["seq"])
        pep_index = peptides.select(
            ["peptide_id", "spectrum_id", "n_reads", "alc", "alc_rank"]
        ).with_row_index("pep_row")

        keyed = (
            scored.join(self._clone_index, on="clone_id", how="left")
            .join(pep_index, on="peptide_id", how="left")
        )
        gamma = (
            self._responsibility(keyed) if competition
            else np.ones(keyed.height)
        )
        confirmed, contradiction = self._cover(
            keyed, gamma, pep_code, pep_offset, pep_length,
            aggregation=aggregation if competition else "alignment_union",
        )
        return self._assemble(confirmed, contradiction, use_prior)

    # -- competition -----------------------------------------------------
    @staticmethod
    def _responsibility(keyed: pl.DataFrame) -> np.ndarray:
        """Posterior that each alignment is the true explanation of its spectrum.

        A softmax over every (reading, clonotype, offset) the spectrum could
        have come from, against a null component of weight 1. An alignment must
        beat the null on its own evidence to earn responsibility, so a weak
        promiscuous hit gets almost none rather than having to divide the
        spectrum's mass with its rivals.

        ``alc`` is **NaN**, not null, on ``db`` rows -- 136 peptides in
        ``dataset_mixa`` and 592 in ``dataset_single_mabs``, though none in
        ``dataset_mixa_gel``, which is why this was invisible on the first
        dataset tried. A single NaN propagates through ``max().over()`` to its
        whole spectrum group and from there, via the coverage union, to most of
        the repertoire: it silently NaN-ed 31,963 of 33,052 clonotypes and sent
        all 12 targets to the bottom of the ranking. NaN is therefore mapped to
        missing, and the result is asserted finite rather than trusted.
        """
        log_weight = (
            pl.col("llr")
            + 2.0 * pl.col("alc").fill_nan(None).fill_null(_ALC_DEFAULT)
            .clip(1e-3, 1.0).log()
            - pl.col("alc_rank").cast(pl.Float64).clip(1.0, None).log()
            - pl.col("n_reads").cast(pl.Float64).clip(1.0, None).log()
        )
        peak = log_weight.max().over("spectrum_id")
        shifted = (log_weight - peak).exp()
        # The null carries weight 1 before shifting, hence exp(-peak) after.
        gamma = keyed.select(
            gamma=shifted / (shifted.sum().over("spectrum_id") + (-peak).exp())
        )["gamma"].to_numpy()
        if not np.isfinite(gamma).all():
            raise ValueError(
                f"{int((~np.isfinite(gamma)).sum())} of {gamma.size} alignment "
                "responsibilities are not finite; check alc/alc_rank/n_reads "
                "for NaN before scoring"
            )
        return gamma

    # -- coverage --------------------------------------------------------
    def _cover(self, keyed: pl.DataFrame, gamma: np.ndarray, pep_code,
               pep_offset, pep_length, aggregation: str = "alignment_union"):
        """Per-residue confirmation and per-clonotype contradiction."""
        n_slots = self.background.code.size
        log_miss = np.zeros(n_slots)
        contradiction = np.zeros(len(self._clone_index))

        clone_row_all = keyed["clone_row"].to_numpy()
        llr_all = keyed["llr"].to_numpy()
        negative = np.minimum(llr_all, 0.0) * gamma
        exclusive = aggregation == "exclusive_spectrum"
        if exclusive:
            # Keep one spectrum contiguous, but still bound residue expansion
            # by _CHUNK. A boundary spectrum's sparse partial sums survive to
            # the next chunk, including a spectrum larger than _CHUNK itself.
            order = np.argsort(keyed["spectrum_id"].to_numpy(), kind="stable")
            keyed = keyed[order]
            gamma = gamma[order]
            spectrum = keyed["spectrum_id"].rle_id().to_numpy()
            carry_spectrum = np.empty(0, dtype=np.int64)
            carry_slot = np.empty(0, dtype=np.int64)
            carry_mass = np.empty(0)
        for begin in range(0, keyed.height, _CHUNK):
            chunk = keyed.slice(begin, _CHUNK)
            size = chunk.height
            pep_row = chunk["pep_row"].to_numpy().astype(np.int64)
            clone_row = chunk["clone_row"].to_numpy()
            start = chunk["start"].to_numpy().astype(np.int64)
            g = gamma[begin:begin + size]

            alignment_of, pep_slot, ref_slot = expand_alignments(
                pep_offset, pep_length, pep_row,
                self.background.offset, clone_row, start,
            )
            agree = (
                (pep_code[pep_slot] == self.background.code[ref_slot])
                & (pep_code[pep_slot] >= 0)
            )
            if exclusive:
                spectra = np.concatenate((carry_spectrum,
                    spectrum[begin + alignment_of[agree]]))
                slots = np.concatenate((carry_slot, ref_slot[agree]))
                mass = np.concatenate((carry_mass, g[alignment_of[agree]]))
                order = np.lexsort((slots, spectra))
                spectra, slots, mass = spectra[order], slots[order], mass[order]
                if mass.size:
                    starts = np.flatnonzero(np.r_[True,
                        (spectra[1:] != spectra[:-1]) | (slots[1:] != slots[:-1])])
                    mass = np.add.reduceat(mass, starts)
                    spectra, slots = spectra[starts], slots[starts]
                complete = spectra != spectrum[begin + size - 1]
                if begin + size == keyed.height:
                    complete = np.ones(mass.size, dtype=bool)
                np.add.at(log_miss, slots[complete],
                          np.log(np.maximum(1.0 - mass[complete], 1e-12)))
                carry_spectrum = spectra[~complete]
                carry_slot = slots[~complete]
                carry_mass = mass[~complete]
            else:
                contribution = np.log(np.maximum(1.0 - g, 1e-12))[alignment_of]
                np.add.at(log_miss, ref_slot[agree], contribution[agree])

        # Contradictions are per clonotype, averaged over the alignment's
        # length so a long peptide is not counted once per residue.
        np.add.at(contradiction, clone_row_all, negative)
        return 1.0 - np.exp(log_miss), np.maximum(contradiction, -_D_MAX)

    # -- assembly --------------------------------------------------------
    def _assemble(self, confirmed, contradiction, use_prior: bool) -> pl.DataFrame:
        bounds = np.append(
            self.background.offset, self.background.code.size
        ).astype(np.int64)
        per_slot = confirmed * self._information
        gained = np.add.reduceat(per_slot, bounds[:-1])
        junctional = np.add.reduceat(
            np.where(self._is_junction, per_slot, 0.0), bounds[:-1]
        )
        out = self.clonotypes.select("clone_id").with_columns(
            confirmed_nats=pl.Series(gained),
            junction_nats=pl.Series(junctional),
            contradiction=pl.Series(contradiction),
            evidence=pl.Series(gained + contradiction),
            # Material confirmed evidence, none of it junctional: the score is
            # about the clonotype's V gene, not the clonotype. Reported, never
            # dropped -- whether such a row is usable is the caller's decision.
            framework_only=pl.Series(
                (gained >= _FLAG_MIN_NATS) & (junctional == 0.0)
            ),
        )
        if not use_prior:
            return _rank_order(out.with_columns(score=pl.col("evidence")))
        isotypes = self.clonotypes.select(isotype_of(self.clonotypes))["isotype"].to_numpy()
        return _rank_order(out.with_columns(
            score=pl.col("evidence") + pl.Series(isotype_log_prior(isotypes))
        ))


def compete_identity_classes(scored: pl.DataFrame, peptides: pl.DataFrame,
                             identities: pl.DataFrame) -> pl.DataFrame:
    """Experimental source-class competition; does not change default scoring.

    ``identities`` maps every clone_id to an explicitly defined identity_class.
    One best placement per spectrum/class prevents sequence aliases and offsets
    from multiplying evidence. Scores compare that placement with all other
    classes and the existing unit-weight null. Exact winning ties stay unresolved.
    This benchmark endpoint is not a calibrated probability of whole-antibody
    identity; downstream operating points require independent calibration.
    """
    if identities['clone_id'].n_unique() != identities.height:
        raise ValueError('identity mapping must have unique clone IDs')
    if identities['identity_class'].null_count():
        raise ValueError('identity classes must be explicit')
    if scored.join(identities.select('clone_id'), on='clone_id', how='anti').height:
        raise ValueError('scored candidate lacks identity class')
    metadata = spectrum_key(peptides).select(
        'peptide_id', 'spectrum_id', 'alc', 'alc_rank', 'n_reads')
    if metadata['peptide_id'].n_unique() != metadata.height:
        raise ValueError('peptide observation IDs must be unique')
    if scored.join(metadata.select('peptide_id'), on='peptide_id', how='anti').height:
        raise ValueError('scored peptide lacks observation metadata')
    keyed = scored.join(identities.select('clone_id', 'identity_class'),
                        on='clone_id', validate='m:1').join(
        metadata, on='peptide_id', validate='m:1')
    log_weight = (pl.col('llr')
        + 2 * pl.col('alc').fill_nan(None).fill_null(_ALC_DEFAULT).clip(1e-3, 1.).log()
        - pl.col('alc_rank').cast(pl.Float64).clip(1., None).log()
        - pl.col('n_reads').cast(pl.Float64).clip(1., None).log())
    result = keyed.with_columns(log_weight.alias('log_weight')).group_by(
        ['spectrum_id', 'identity_class']).agg(pl.col('log_weight').max())
    if result.filter(~pl.col('log_weight').is_finite() | pl.col('log_weight').is_null()).height:
        raise ValueError('nonfinite source-class evidence')
    result = result.with_columns(
        pl.col('log_weight').max().over('spectrum_id').clip(lower_bound=0.).alias('_peak'))
    result = result.with_columns((pl.col('log_weight')-pl.col('_peak')).exp().alias('_weight'))
    result = result.with_columns(
        (pl.col('_weight') / (pl.col('_weight').sum().over('spectrum_id') +
          (-pl.col('_peak')).exp())).alias('responsibility'),
        (pl.col('log_weight') == pl.col('log_weight').max().over('spectrum_id')).alias('_best'))
    return result.with_columns(
        (pl.col('_best') & (pl.col('_best').sum().over('spectrum_id') == 1) &
         (pl.col('log_weight') > 0)).alias('unique_map_above_null')
    ).drop('_peak', '_weight', '_best').sort(['spectrum_id', 'identity_class'])


def fit_identity_mixture(competition: pl.DataFrame, peptides: pl.DataFrame, *,
                         max_iterations: int = 1000, tolerance: float = 1e-8):
    """Experimental sparse conditional mixture with explicit unit-likelihood null.

    Fit nonnegative class weights from observed likelihoods only. Repeated I/L
    peptide species share total unit weight, preserving distinct spectra while
    preventing reacquisition from multiplying fitting evidence. Multiple sources
    can receive positive weights; per-spectrum MAP selection is never required.
    Removal contrasts keep other fitted weights fixed and renormalize them; they
    are not profiled likelihood ratios or calibrated source-presence probabilities.
    Only matched classes are returned; callers must left join their full declared
    identity catalog and assign zero to genuinely unmatched classes. Exact equal
    sparse likelihood profiles form observational equivalence groups; individual
    members are not identifiable and must not become unique source calls. Unavailable
    contrasts and nonconverged fixed-other results alone never establish presence.
    The separate certified-profile API can bound presence evidence without
    assuming that the full or reduced fitting point is an optimum.
    """
    if type(max_iterations) is not int or max_iterations < 1:
        raise ValueError('max_iterations must be a positive integer')
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError('tolerance must be finite and positive')
    if competition.select(pl.struct('spectrum_id','identity_class').n_unique()).item() != competition.height:
        raise ValueError('competition rows must be unique')
    if competition.filter(~pl.col('log_weight').is_finite() | pl.col('log_weight').is_null()).height:
        raise ValueError('class likelihoods must be finite')
    if peptides['seq'].null_count():
        raise ValueError('mixture peptide sequences must be explicit')
    if competition['spectrum_id'].null_count() or competition['identity_class'].null_count():
        raise ValueError('competition spectrum and class identities must be explicit')
    species=spectrum_key(peptides).select('spectrum_id',
        pl.col('seq').str.replace_all('I','L',literal=True).alias('species')).unique()
    events=species.group_by('spectrum_id').agg(pl.col('species').sort().alias('species'))
    # A spectrum may have several alternative peptide calls; their whole sorted
    # set defines its observation species rather than multiplying the event.
    events=events.with_columns(pl.col('species').list.join('|').alias('species_key')).sort('spectrum_id')
    events=events.with_columns((1./pl.len().over('species_key')).alias('weight')).with_row_index('event_row')
    # Missing edges mean zero likelihood. Finite explicit log likelihoods are
    # positive, so sorted sparse profiles exactly characterize equality without
    # materializing an observations-by-classes matrix or inventing a tolerance.
    classes=competition.group_by('identity_class').agg(
        pl.struct('spectrum_id','log_weight').sort_by('spectrum_id').alias('_profile')
    ).sort('identity_class').with_row_index('class_row')
    classes=classes.with_columns(
        pl.col('class_row').min().over('_profile').alias('identity_equivalence_group'),
        pl.len().over('_profile').alias('identity_equivalence_size')
    ).drop('_profile').with_columns(
        (pl.col('identity_equivalence_size')==1).alias('identifiable'))
    edges=competition.join(events.select('spectrum_id','event_row'),on='spectrum_id',how='left',validate='m:1').join(classes,on='identity_class',validate='m:1').sort(['event_row','class_row'])
    if edges['event_row'].null_count():raise ValueError('class likelihood lacks observed spectrum')
    er=edges['event_row'].to_numpy();cr=edges['class_row'].to_numpy()
    ne,nc=events.height,classes.height
    ew=events['weight'].to_numpy();total=ew.sum()
    if not ne or not total:raise ValueError('mixture requires observed spectra')
    logl=edges['log_weight'].to_numpy();shift=np.zeros(ne)
    np.maximum.at(shift,er,logl)
    likelihood=np.exp(logl-shift[er]);null=np.exp(-shift)
    weights=np.full(nc+1,1./(nc+1));history=[];converged=False
    log_null_weight=-np.log(nc+1.)
    for iteration in range(max_iterations):
        denominator=weights[-1]*null+np.bincount(er,weights=weights[cr]*likelihood,minlength=ne)
        if (denominator<=0).any():raise ValueError('zero mixture likelihood')
        resp=weights[cr]*likelihood/denominator[er]
        updated=np.r_[np.bincount(cr,weights=ew[er]*resp,minlength=nc),
            np.sum(ew*weights[-1]*null/denominator)]/total
        log_null_weight += float(np.logaddexp.reduce(np.log(ew)-shift-np.log(denominator))-np.log(total))
        objective=float(np.sum(ew*(np.log(denominator)+shift)))
        if not np.isfinite(objective):
            raise ValueError('nonfinite mixture objective')
        history.append(objective)
        delta=float(np.max(np.abs(updated-weights)));weights=updated
        if delta<=tolerance:converged=True;break
    denominator=weights[-1]*null+np.bincount(er,weights=weights[cr]*likelihood,minlength=ne)
    final_objective=float(np.sum(ew*(np.log(denominator)+shift)))
    gradient=np.r_[np.bincount(cr,weights=ew[er]*likelihood/denominator[er],minlength=nc),
        np.sum(ew*null/denominator)]
    if not np.isfinite(gradient).all() or not np.isfinite(final_objective):
        raise ValueError('nonfinite mixture gradient/objective; no certified bound')
    max_gradient=float(gradient.max())
    dual_gap=max(0.,max_gradient-float(total))
    # Work in log space for leave-one-class-out denominators: subtracting a
    # dominant contribution from the full denominator loses every other source.
    with np.errstate(divide='ignore'):
        edge_log = np.log(weights[cr]) + logl
        null_log = np.full(ne, log_null_weight)
    all_event = np.r_[er, np.arange(ne)]
    all_log = np.r_[edge_log, null_log]
    peak = np.full(ne, -np.inf)
    np.maximum.at(peak, all_event, all_log)
    indices = np.arange(all_log.size)
    first = np.full(ne, all_log.size, dtype=np.int64)
    np.minimum.at(first, all_event, np.where(all_log == peak[all_event], indices, all_log.size))
    remaining = indices != first[all_event]
    second = np.full(ne, -np.inf)
    np.maximum.at(second, all_event[remaining], all_log[remaining])
    finite = remaining & np.isfinite(all_log)
    rest_sum = np.bincount(all_event[finite],
        weights=np.exp(all_log[finite] - second[all_event[finite]]), minlength=ne)
    with np.errstate(divide='ignore'):
        rest_log = second + np.log(rest_sum)
    full_log = np.logaddexp(peak, rest_log)
    # For every nonmaximum term its share is at most 1/2, so subtraction here
    # is well conditioned; only the maximum needs the separately summed rest.
    with np.errstate(divide='ignore', invalid='ignore'):
        leave_log = full_log[er] + np.log1p(-np.exp(edge_log-full_log[er]))
    winner = indices[:er.size] == first[er]
    leave_log[winner] = rest_log[er[winner]]
    change = full_log[er] - leave_log
    change[np.isneginf(edge_log)] = 0.
    contrast = np.bincount(cr, weights=ew[er]*change, minlength=nc).astype(float)
    with np.errstate(divide='ignore', invalid='ignore'):
        log_other = np.log1p(-weights[:-1])
        component_logs = np.r_[np.log(weights[:-1]), log_null_weight]
    dominant = int(np.argmax(component_logs))
    if dominant < nc:
        log_other[dominant] = np.logaddexp.reduce(np.delete(component_logs, dominant))
    contrast += total*log_other
    available = np.isfinite(contrast)
    contrast[~available] = np.nan
    result=classes.drop('class_row').with_columns(pl.Series('mixture_weight',weights[:-1]),
        pl.Series('fixed_other_removal_contrast',contrast),
        pl.Series('removal_contrast_available',available))
    return result,dict(converged=converged,iterations=iteration+1,null_weight=float(weights[-1]),
        objective=history,final_objective=final_objective,dual_gap=dual_gap,
        max_component_gradient=max_gradient,weighted_observations=float(total),observations=ne,
        classes=nc,placement_edges=edges.height,max_weight_change=delta,
        unavailable_contrasts=int((~available).sum()),log_null_weight=float(log_null_weight),
        observational_equivalence_groups=classes['identity_equivalence_group'].n_unique(),
        individually_unidentifiable_classes=classes.filter(~pl.col('identifiable')).height,
        rule='conditional mixture; fixed-other renormalized removal contrast, not profiled LR')


def profile_identity_removals(competition: pl.DataFrame, peptides: pl.DataFrame, *,
                              minimum_contrast: float = 0.,
                              classes_to_profile=None,
                              max_iterations: int = 1000, tolerance: float = 1e-8):
    """Conservative experimental presence evidence from refitted class removals.

    For concave weighted log likelihood Q on the weight simplex, the gradient
    gap U=max(gradient)-sum(observation weights) bounds optimization error. The
    reported lower bound max(0,Q(full fitted)-Q(reduced fitted)-U(reduced))
    cannot exceed the exact profiled removal loss. No likelihood-ratio null
    distribution or posterior calibration is implied. Full-fit optimization
    error is included in the true profiled-loss upper bound. The conservative
    score itself is additionally bounded above by max(0,fixed-other contrast),
    since the reduced optimum dominates the feasible fixed-other removal. Either
    valid upper bound can screen refits below the requested operating threshold.

    ``classes_to_profile`` explicitly bounds development/profiling calls; omitted
    classes remain unavailable, never assigned an accepted zero/error surrogate.
    Equal observed profiles remain individually unresolved regardless of score.
    """
    if not np.isfinite(minimum_contrast) or minimum_contrast < 0:
        raise ValueError('minimum_contrast must be finite and nonnegative')
    full,receipt=fit_identity_mixture(competition,peptides,
        max_iterations=max_iterations,tolerance=tolerance)
    requested=set(full['identity_class']) if classes_to_profile is None else set(classes_to_profile)
    unknown=requested-set(full['identity_class'])
    if unknown:raise ValueError('requested profile class lacks a fitted likelihood')
    rows=[];refits=[];failures=[]
    for row in full.iter_rows(named=True):
        identity=row['identity_class']
        fixed=row['fixed_other_removal_contrast']
        upper=max(0.,fixed+receipt['dual_gap']) if row['removal_contrast_available'] else float('nan')
        score_upper=max(0.,fixed) if row['removal_contrast_available'] else float('nan')
        lower=None;reduced_objective=None;reduced_gap=None
        if not row['identifiable']:
            status='observationally_equivalent'
        elif not row['removal_contrast_available']:
            status='fixed_other_unavailable'
        elif upper<=minimum_contrast or score_upper<=minimum_contrast:
            lower=0.;status='screened_upper_bound'
        elif identity not in requested:
            status='not_profiled'
        else:
            try:
                _,reduced_receipt=fit_identity_mixture(
                    competition.filter(pl.col('identity_class')!=identity),peptides,
                    max_iterations=max_iterations,tolerance=tolerance)
            except ValueError as error:
                status='reduced_fit_error'
                failures.append(dict(identity_class=identity,error=str(error)))
            else:
                reduced_objective=reduced_receipt['final_objective']
                reduced_gap=reduced_receipt['dual_gap']
                lower=max(0.,receipt['final_objective']-reduced_objective-reduced_gap)
                if np.isfinite(lower):
                    status='profiled_certified_bound'
                else:
                    lower=None;status='reduced_bound_error'
                    failures.append(dict(identity_class=identity,error='nonfinite certified removal bound'))
                refits.append(dict(identity_class=identity,**reduced_receipt))
        rows.append(dict(identity_class=identity,profile_status=status,
            profiled_removal_lower_bound=lower,profiled_removal_upper_bound=upper,
            certified_score_upper_bound=score_upper,
            reduced_final_objective=reduced_objective,reduced_dual_gap=reduced_gap,
            profile_score_available=lower is not None))
    bounds=pl.DataFrame(rows,schema_overrides={
        'identity_class':full.schema['identity_class'],
        'profiled_removal_lower_bound':pl.Float64,'profiled_removal_upper_bound':pl.Float64,
        'certified_score_upper_bound':pl.Float64,
        'reduced_final_objective':pl.Float64,'reduced_dual_gap':pl.Float64,
        'profile_status':pl.String,'profile_score_available':pl.Boolean}) if rows else full.select('identity_class').with_columns(
            pl.lit(None,dtype=pl.String).alias('profile_status'),
            *[pl.lit(None,dtype=pl.Float64).alias(name) for name in
              ['profiled_removal_lower_bound','profiled_removal_upper_bound','certified_score_upper_bound',
               'reduced_final_objective','reduced_dual_gap']],
            pl.lit(False).alias('profile_score_available'))
    return full.join(bounds,on='identity_class',validate='1:1'),dict(
        full_fit=receipt,reduced_fits=refits,reduced_failures=failures,minimum_contrast=float(minimum_contrast),
        requested_profiles=len(requested),performed_profiles=len(refits),
        rule='max(0,Q(full fitted)-Q(reduced fitted)-dual_gap(reduced)); conservative profiled-removal lower bound')


def compete_identity_emissions(components: pl.DataFrame, peptides: pl.DataFrame,
                               identities: pl.DataFrame):
    """Optional common-null class competition from observed-AA log emissions.

    The null is the uniform average of each declared class's best emission for
    that peptide, with absent search-supported placements assigned zero. This
    truncated candidate-pool background is not a calibrated false-PSM density.
    Original candidate-context LLRs are not used as comparable source densities.
    Return competition plus per-peptide null construction/provenance metadata.
    """
    if identities['clone_id'].n_unique()!=identities.height or identities['identity_class'].null_count():
        raise ValueError('identity mapping must be unique and explicit')
    if components.join(identities.select('clone_id'),on='clone_id',how='anti').height:
        raise ValueError('scored candidate lacks identity class')
    if components.select('peptide_id').unique().join(
            peptides.select('peptide_id'),on='peptide_id',how='anti').height:
        raise ValueError('scored peptide lacks observation metadata')
    n_classes=identities['identity_class'].n_unique()
    if not n_classes:raise ValueError('common null requires a full declared class catalog')
    value=pl.col('log_emission')
    if components.filter(value.is_null() | value.is_nan() | (value==float('inf'))).height:
        raise ValueError('observed-AA emission unavailable/nonfinite')
    # -inf is a genuine zero emission (for example a mismatch with zero error).
    supported=components.filter(value.is_finite())
    class_emissions=supported.join(identities.select('clone_id','identity_class'),on='clone_id',validate='m:1').group_by(
        ['peptide_id','identity_class']).agg(pl.col('log_emission').max())
    class_emissions=class_emissions.with_columns(
        pl.col('log_emission').max().over('peptide_id').alias('_peak'))
    common=class_emissions.group_by('peptide_id').agg(
        pl.col('_peak').first(),
        (pl.col('log_emission')-pl.col('_peak')).exp().sum().alias('_sum'),
        pl.len().alias('supporting_classes')).with_columns(
            (pl.col('_peak')+pl.col('_sum').log()-np.log(n_classes)).alias('common_null_log_emission')
        ).drop('_peak','_sum')
    null=peptides.select('peptide_id').join(common,on='peptide_id',how='left',validate='1:1').with_columns(
        pl.col('supporting_classes').fill_null(0),pl.lit(n_classes).alias('declared_classes')
    ).with_columns(pl.when(pl.col('supporting_classes')>0).then(pl.lit('common_class_average')).otherwise(
        pl.lit('null_only_no_supported_class')).alias('null_status')).sort('peptide_id')
    relative=supported.join(common.select('peptide_id','common_null_log_emission'),on='peptide_id',validate='m:1').with_columns(
        (pl.col('log_emission')-pl.col('common_null_log_emission')).alias('llr'))
    return compete_identity_classes(relative,peptides,identities),null


_RNA_SUPPORT_UNITS = 'qualifying native raw/assembled observations; not independent molecules'


def identity_rna_support(candidates: pl.DataFrame, identities: pl.DataFrame) -> pl.DataFrame:
    """Full-class maximum native-root observation support, without alias sums.

    Parent counts concern qualifying native observations, which can include raw
    and assembled duplicates. They are not independent molecular counts, full-AA
    sequence quality, protein amounts or physical-source labels. Zero/fallback
    support remains in the full declared class population.
    """
    required={'clone_id','parent_clone_id','cdr3_seq','representative_policy','native_root_matching_observations'}
    if not required<=set(candidates.columns) or not {'clone_id','identity_class'}<=set(identities.columns):
        raise ValueError('native-root candidate and full identity metadata are required')
    if candidates['representative_policy'].null_count() or set(candidates['representative_policy'])!={'native_root'}:
        raise ValueError('RNA reliability requires the native_root representative policy')
    if candidates['clone_id'].null_count() or candidates['clone_id'].n_unique()!=candidates.height or identities['clone_id'].null_count() or identities['clone_id'].n_unique()!=identities.height:
        raise ValueError('candidate and identity clone IDs must be unique and explicit')
    if candidates.select('clone_id').join(identities.select('clone_id'),on='clone_id',how='anti').height or identities.select('clone_id').join(candidates.select('clone_id'),on='clone_id',how='anti').height:
        raise ValueError('identity mapping must cover the complete candidate population')
    count='native_root_matching_observations'
    if not candidates.schema[count].is_integer() or candidates[count].null_count() or candidates.filter((pl.col(count)<0)|(pl.col(count)>2**53)).height:
        raise ValueError('RNA support must be nonmissing nonnegative native integer observation counts')
    if candidates['parent_clone_id'].null_count() or candidates['cdr3_seq'].null_count() or identities['identity_class'].null_count():
        raise ValueError('native parents and junction identities must be explicit')
    keyed=candidates.select('clone_id','parent_clone_id','cdr3_seq',count).join(identities.select('clone_id','identity_class'),on='clone_id',validate='1:1')
    if keyed.filter(pl.col('identity_class')!=pl.col('cdr3_seq').str.replace_all('I','L',literal=True)).height:
        raise ValueError('identity mapping differs from the observed I/L-equivalent junction')
    parent=keyed.group_by('parent_clone_id').agg(pl.col('identity_class').n_unique().alias('_classes'),
        pl.col(count).n_unique().alias('_counts'))
    if parent.filter((pl.col('_classes')!=1)|(pl.col('_counts')!=1)).height:
        raise ValueError('aliases of one native parent must agree on junction class and observation count')
    return keyed.group_by(['parent_clone_id','identity_class']).agg(pl.col(count).max()).group_by('identity_class').agg(
        pl.col(count).max().cast(pl.Int64).alias('rna_root_observations'),pl.len().alias('native_parents')
    ).with_columns(pl.lit(_RNA_SUPPORT_UNITS).alias('rna_support_units')).sort('identity_class')


def apply_rna_junction_reliability(class_scores: pl.DataFrame, support: pl.DataFrame,
                                   knots, q) -> pl.DataFrame:
    """Append S=max(0,B+log(q)) from frozen RNA calibration, with no truth input.

    Interpolation is linear between unique support knots, constant outside them.
    All original rows/columns, ambiguity and unavailable states are retained.
    q=1 permits exact baseline parity; fitted smoothed curves lie strictly below
    one. This only penalizes a conservative score and is not a protein posterior
    or FDR estimate. Missing RNA metadata fails closed rather than dropping rows.
    """
    x=np.asarray(knots,dtype=np.float64);prob=np.asarray(q,dtype=np.float64)
    if x.ndim!=1 or prob.ndim!=1 or not x.size or x.size!=prob.size or not np.isfinite(x).all() or not np.isfinite(prob).all() or (x<0).any() or (np.diff(x)<=0).any() or ((prob<=0)|(prob>1)).any() or (np.diff(prob)<0).any():
        raise ValueError('frozen RNA curve requires increasing finite support and monotone probabilities in (0,1]')
    required={'identity_class','rna_root_observations','rna_support_units'}
    if not required<=set(support.columns):raise ValueError('full RNA support metadata and observation units are required')
    for frame in (class_scores,support):
        if frame['identity_class'].null_count() or frame['identity_class'].n_unique()!=frame.height:
            raise ValueError('class tables must have unique explicit identities')
    if class_scores.select('identity_class').join(support.select('identity_class'),on='identity_class',how='anti').height or support.select('identity_class').join(class_scores.select('identity_class'),on='identity_class',how='anti').height:
        raise ValueError('RNA support must preserve the full score-class population')
    if support['rna_support_units'].null_count() or support.filter(pl.col('rna_support_units')!=_RNA_SUPPORT_UNITS).height:
        raise ValueError('RNA support units must be qualifying native observations, not molecular abundance')
    if not support.schema['rna_root_observations'].is_integer() or support['rna_root_observations'].null_count() or support.filter((pl.col('rna_root_observations')<0)|(pl.col('rna_root_observations')>2**53)).height:
        raise ValueError('RNA support counts must be explicit nonnegative integers')
    if class_scores['profile_score_available'].null_count():raise ValueError('profile availability must be explicit')
    value=pl.col('profiled_removal_lower_bound')
    if class_scores.filter(pl.col('profile_score_available') & (value.is_null()|~value.is_finite()|(value<0))).height:
        raise ValueError('available profile scores must be finite nonnegative bounds')
    out=class_scores.with_row_index('_rna_row').join(support,on='identity_class',validate='1:1').sort('_rna_row').drop('_rna_row')
    reliability=np.interp(out['rna_root_observations'].to_numpy(),x,prob)
    return out.with_columns(pl.Series('rna_junction_reliability',reliability)).with_columns(
        pl.when(pl.col('profile_score_available')).then((value+pl.col('rna_junction_reliability').log()).clip(lower_bound=0.)).otherwise(None).alias('rna_adjusted_removal_lower_bound'),
        pl.col('profile_score_available').alias('rna_adjusted_score_available'))
