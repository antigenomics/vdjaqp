"""Simulate library sampling from an explicitly fitted empirical repertoire.

This conditional multinomial model matches observed molecule depth and junction
frequencies. It is fitted mimicry, not independent RNA/MS pipeline validation;
unseen identities cannot be inferred. All replicate draws are generated afresh.
"""
import argparse
import hashlib
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np
import polars as pl


def draw_library(task):
    counts, depth, replicates, seed = task
    probabilities = np.r_[counts, depth-counts.sum()] / depth
    assert np.all(probabilities >= 0) and np.isclose(probabilities.sum(), 1)
    return np.random.default_rng(seed).multinomial(depth, probabilities, size=replicates)[:, :-1].astype(np.uint32)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--abundance', required=True, type=Path)
    p.add_argument('--processing', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--replicates', type=int, default=200)
    p.add_argument('--seed', type=int, default=20261116)
    a = p.parse_args()
    if a.replicates < 100:
        p.error('at least 100 replicates are required for the displayed interval')
    a.output.mkdir(parents=True, exist_ok=False)
    source = pl.read_parquet(a.abundance)
    processing = pl.read_csv(a.processing, separator='\t').filter(pl.col('input')!='joint')
    matrix = source.pivot(on='subset_id', index='corrected_junction_nt_key', values='molecule_groups').fill_null(0).sort('corrected_junction_nt_key')
    libraries = processing['input'].to_list()
    counts = np.column_stack([matrix[s].to_numpy() for s in libraries])
    depths = processing['molecule_groups'].to_numpy()
    # Three coarse independent native draws; no per-identity jobs or nested pool.
    with mp.get_context('spawn').Pool(min(len(libraries), mp.cpu_count())) as pool:
        draws = pool.map(draw_library, [(counts[:,i], int(depths[i]), a.replicates, a.seed+i) for i in range(len(libraries))])
    curves, summaries = [], []
    for i, (library, simulated) in enumerate(zip(libraries, draws)):
        observed = np.sort(counts[:,i])[::-1]
        rank = np.unique(np.r_[1, np.geomspace(1, len(observed), 180).astype(int), (observed>0).sum()])
        cumulative = np.cumsum(np.sort(simulated,axis=1)[:,::-1],axis=1)[:,rank-1] / depths[i]
        lower, median, upper = np.quantile(cumulative,[.025,.5,.975],axis=0)
        for j,r in enumerate(rank):
            curves.append(dict(library=library,rank=int(r),observed=float(observed[:r].sum()/depths[i]),
                               simulated_median=median[j],lower=lower[j],upper=upper[j]))
        richness = (simulated>0).sum(axis=1)
        summaries.append(dict(quantity='detected_junctions',library=library,observed=int((counts[:,i]>0).sum()),
                              median=float(np.median(richness)),lower=float(np.quantile(richness,.025)),upper=float(np.quantile(richness,.975))))
    membership = sum((x>0).astype(np.uint8) for x in draws)
    weights = sum(x.astype(np.uint64) for x in draws)
    observed_shared = (counts>0).sum(axis=1)>1
    for unit, values, obs in [
        ('shared_identity_fraction',(membership>1).sum(axis=1)/(membership>0).sum(axis=1),observed_shared.mean()),
        ('shared_molecule_fraction',np.where(membership>1,weights,0).sum(axis=1)/weights.sum(axis=1),counts[observed_shared].sum()/counts.sum())]:
        summaries.append(dict(quantity=unit,library='joint',observed=float(obs),median=float(np.median(values)),
                              lower=float(np.quantile(values,.025)),upper=float(np.quantile(values,.975))))
    pl.DataFrame(curves).write_csv(a.output/'curves.tsv',separator='\t')
    pl.DataFrame(summaries).write_csv(a.output/'summary.tsv',separator='\t')
    (a.output/'receipt.json').write_text(json.dumps({
        'model':'empirical per-library multinomial with unresolved remainder category',
        'scope':'conditional sampling mimicry; no unseen identities, biological replication, read errors, assembly or MS simulation',
        'replicates':a.replicates,'seed':a.seed,'molecule_depths':dict(zip(libraries,map(int,depths))),
        'junctions':matrix.height,'interval':'pointwise central 95% simulation interval, not a biological confidence interval',
        'input_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [a.abundance,a.processing]},
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    },indent=2)+'\n')


if __name__ == '__main__':
    main()
