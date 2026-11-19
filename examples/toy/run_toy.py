"""Purpose: a self-contained, fully synthetic walk-through of VDJ-AQP scoring.

Context: tutorial and smoke example only. It builds a small invented heavy-chain
repertoire (one target, one junction relative, one CDR2 variant and 300 background
sequences sharing the same IGHV3-23 framework) and five tryptic peptides from the
target, then runs the unchanged search, background, error-model and aggregation path
(`Evidence.from_frames`). No private data, calibration or scientific result is
involved; the fixed error model is the documented database-search assumption.

Usage: .venv/bin/python examples/toy/run_toy.py [--output DIR]
"""
import argparse
from pathlib import Path

import numpy as np
import polars as pl

from vdjaqp.errormodel import ErrorModel, mass_prior
from vdjaqp.pipeline import Evidence

# Supplied region labels: 0 FR1, 1 CDR1, 2 FR2, 3 CDR2, 4 FR3, 5 CDR3 anchors, 6 junction interior, 7 FR4.
FR1, CDR1, FR2 = 'EVQLLESGGGLVQPGGSLRLSCAAS', 'GFTFSSYA', 'MSWVRQAPGKGLEWVS'
CDR2, FR3, FR4 = 'AISGSGGST', 'YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYY', 'GQGTLVTVSS'
AMINO = list('ACDEFGHKLMNPQRSTVWY')
PEPTIDES = ('LSCAASGFTFSSYAMSWVR',      # FR1-CDR1-FR2: shared by every candidate
            'GLEWVSAISGSGGSTYYADSVK',   # FR2-CDR2-FR3: shared except the CDR2 variant
            'NTLYLQMNSLR',              # FR3: shared by every candidate
            'AEDTAVYYCAKDR',            # FR3 into the junction start
            'GYSSGWYFDYWGQGTLVTVSS')    # junction interior into FR4: distinguishing


def record(clone_id, junction, cdr2=CDR2):
    """One candidate with per-residue region labels; the junction includes C and W anchors."""
    parts = [(FR1, 0), (CDR1, 1), (FR2, 2), (cdr2, 3), (FR3, 4), (junction, 5), (FR4, 7)]
    seq = ''.join(p for p, _ in parts)
    regions = [r for p, r in parts for _ in p]
    start = len(FR1 + CDR1 + FR2 + cdr2 + FR3)
    for i in range(start + 3, start + len(junction) - 3):
        regions[i] = 6
    return dict(clone_id=clone_id, seq=seq, region_ids=regions, cdr3_start=start,
                cdr3_end=start + len(junction), cdr3_seq=junction,
                v_name='IGHV3-23', j_name='IGHJ4', c_name='IGHG1')


def repertoire(n_background=300, seed=7):
    rng = np.random.default_rng(seed)
    rows = [record('target', 'CAKDRGYSSGWYFDYW'),
            record('junction_relative', 'CAKDRGYSSGWYFDNW'),   # one junction substitution
            record('cdr2_variant', 'CAKEGLTGTTFDYW', cdr2='VISYDGSNK')]
    for i in range(n_background):
        core = ''.join(rng.choice(AMINO, int(rng.integers(6, 12))))
        rows.append(record(f'background_{i:03d}', f'CAK{core}FDYW'))
    return pl.DataFrame(rows)


def peptides():
    n = len(PEPTIDES)
    return pl.DataFrame({
        'peptide_id': [f'toy/{i}/2/50/10/1/{p}' for i, p in enumerate(PEPTIDES)],
        'seq': list(PEPTIDES),
        'local_confidence': [[float('nan')] * len(p) for p in PEPTIDES],
        'search_type': ['db'] * n, 'alc': [float('nan')] * n, 'alc_rank': [1.] * n, 'n_reads': [1] * n})


def run():
    clones, peps = repertoire(), peptides()
    # Fixed database-search assumption: e_db (keyword) is the per-residue error used for
    # database hits, which carry no residue confidence. The fourth positional argument is
    # a substitution count recorded for calibrated models, not an error rate.
    error = ErrorModel(np.array([0., 1.]), np.array([.1, .01]), mass_prior(), 0, e_db=0.10)
    evidence = Evidence.from_frames(clones, peps, error_model=error)
    scores = evidence.score(use_prior=False)
    counts = evidence.match_counts().rename({'score': 'peptide_matches'})
    table = scores.join(counts, on='clone_id').select(
        'clone_id', 'peptide_matches', 'confirmed_nats', 'junction_nats', 'framework_only',
        accepted=pl.col('confirmed_nats') >= 20.0).sort(['confirmed_nats', 'clone_id'], descending=[True, False])
    return clones, peps, evidence, table


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--output', type=Path, help='fresh directory for the score table and inputs')
    args = parser.parse_args()
    clones, peps, evidence, table = run()
    with pl.Config(tbl_rows=8, tbl_width_chars=120, float_precision=3):
        print(f'{clones.height} candidates, {peps.height} peptides, {evidence.matches.height} alignments')
        print(table.head(8))
    print(f"Candidates matched by at least one peptide: {table.filter(pl.col('peptide_matches') > 0).height}; "
          f"accepted at 20 confirmed nats: {table.filter('accepted')['clone_id'].to_list()}")
    if args.output:
        args.output.mkdir(parents=True, exist_ok=False)
        clones.write_parquet(args.output / 'candidates.parquet')
        peps.write_parquet(args.output / 'peptides.parquet')
        table.write_csv(args.output / 'scores.tsv', separator='\t')
