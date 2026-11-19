"""Purpose: score native synthetic PSM exports against observed ARDA candidates.

Context: integration fixture with the existing fixed DB error-rate assumption;
not a calibrated simulation-performance experiment or a historical rerun.
"""
import argparse
import json
from pathlib import Path
import resource
import sys
import time

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidates',type=Path);p.add_argument('psms',type=Path);p.add_argument('output',type=Path);p.add_argument('--run-id',required=True);p.add_argument('--threads',type=int,default=None);a=p.parse_args()
    import vdjaqp
    vdjaqp.resolve_cpu_budget(a.threads)
    import numpy as np
    import polars as pl
    from vdjaqp.pipeline import Evidence
    from vdjaqp.errormodel import ErrorModel, mass_prior, _E_DB_DEFAULT
    from vdjaqp.readout import read_idxml
    start=time.perf_counter();clones=pl.read_parquet(a.candidates)
    peptides=read_idxml(a.psms,run_id=a.run_id)
    model=ErrorModel(np.array([0.,1.]),np.array([.1,.01]),mass_prior(),_E_DB_DEFAULT)
    evidence=Evidence.from_frames(clones,peptides,error_model=model)
    scores=evidence.score()
    assert scores.height==clones.height
    assert np.isfinite(scores['score'].to_numpy()).all()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    scores.write_parquet(a.output,compression='zstd')
    peptides.write_parquet(a.output.with_name(a.output.stem+'_peptides.parquet'),compression='zstd')
    report={'observed_candidates':clones.height,'psm_hits':peptides.height,
        'native_spectra_with_hits':peptides['spectrum_id'].n_unique(),'scored_candidates':scores.height,
        'dropped_candidates':clones.height-scores.height,'alignment_rows':evidence.matches.height,
        'wall_seconds':time.perf_counter()-start,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
        'error_model':'existing fixed DB error-rate assumption; no calibration fitted','e_db':_E_DB_DEFAULT}
    a.output.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
