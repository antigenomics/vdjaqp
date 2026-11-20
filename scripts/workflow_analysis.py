"""Score observed tables with an explicit external model; truth is never loaded."""
import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import time


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('candidates', 'peptides', 'model', 'output'):
        p.add_argument(name, type=Path)
    p.add_argument('--producer', required=True)
    p.add_argument('--cpus', type=int)
    a = p.parse_args()
    import vdjaqp
    cpus = vdjaqp.resolve_cpu_budget(a.cpus)
    import numpy as np
    import polars as pl
    from vdjaqp.adapter import to_airr_compatible, from_airr_compatible
    from vdjaqp.data import to_gene_level
    from vdjaqp.errormodel import ErrorModel
    from vdjaqp.pipeline import Evidence
    from vdjaqp.airr_capture import capture_tables
    start = time.perf_counter()
    if a.output.exists():
        p.error('output must not exist')
    clones = pl.read_parquet(a.candidates)
    if 'background_member' not in clones.columns or clones['background_member'].null_count():
        p.error('explicit background_member is required; never infer from truth')
    clones = from_airr_compatible(to_airr_compatible(clones, source_format=a.producer))
    peptides = pl.read_parquet(a.peptides)
    if 'source_file' in peptides.columns:
        peptides = peptides.drop('source_file')
    model = json.loads(a.model.read_text())
    x, e, substitution = (np.asarray(model[k], dtype=float)
                          for k in ('conf_x', 'conf_e', 'substitution'))
    # Native isotonic calibration retains ties; preserve np.interp tie semantics.
    if (x.ndim != 1 or not x.size or e.shape != x.shape
            or not np.isfinite(x).all() or not np.isfinite(e).all()
            or (np.diff(x) < 0).any() or (x < 0).any() or (x > 1).any()
            or (e <= 0).any() or (e >= 1).any()
            or substitution.shape != (19,19) or not np.isfinite(substitution).all()
            or (substitution < 0).any() or not np.allclose(substitution.sum(axis=1),1)
            or not np.allclose(substitution.diagonal(),0)
            or not 0 < model['e_db'] < 1):
        p.error('invalid confidence/error probabilities or substitution matrix')
    error = ErrorModel(x, e, substitution, model.get('n_substitutions', 0), e_db=model['e_db'])
    evidence = Evidence.from_frames(to_gene_level(clones), peptides, error_model=error)
    scores = evidence.score()
    if scores.height != clones.height or not np.isfinite(scores['score'].to_numpy()).all():
        raise ValueError('incomplete or nonfinite candidate scores')
    a.output.mkdir()
    scores.write_parquet(a.output/'scores.parquet', compression='zstd')
    # No known-answer IDs are inferred for generic input; explicit source roles suffice.
    references=set(clones.filter(~pl.col('background_member'))['clone_id'])
    for name,frame in capture_tables(evidence,scores,references).items():
        frame.write_parquet(a.output/(name+'.parquet'),compression='zstd')
    evidence.match_counts().write_parquet(a.output/'match_counts.parquet', compression='zstd')
    runtime = {name: version(name) for name in ('vdjaqp', 'seqtree', 'numpy', 'polars')}
    runtime['python'] = platform.python_version()
    report = dict(runtime=runtime, producer=a.producer, cpus=cpus, candidates=clones.height,
                  peptides=peptides.height, alignments=evidence.matches.height,
                  seconds=time.perf_counter()-start, truth_loaded=False,
                  inputs={name:dict(sha256=sha256(path))
                          for name,path in [('candidates',a.candidates),('peptides',a.peptides),('model',a.model)]},
                  outputs={path.name: dict(sha256=sha256(path))
                           for path in sorted(a.output.glob('*.parquet'))})
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__ == '__main__':
    main()
