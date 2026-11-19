"""Create a tiny technical acceptance fixture, never a scientific benchmark."""
import json
from pathlib import Path
import sys
import numpy as np
import polars as pl
from vdjaqp.errormodel import mass_prior

out = Path(sys.argv[1]).resolve()
out.mkdir(parents=True, exist_ok=False)
pl.DataFrame({'clone_id':['c0'], 'seq':['ACDEFGHKLMNP'],
    'cdr3_seq':['ACDEFGHKLMNP'], 'region_ids':[[6]*12],
    'cdr3_start':[0], 'cdr3_end':[12], 'background_member':[True],
    'v_name':['IGHV1-2'], 'j_name':['IGHJ4'], 'c_name':['IGHG1']}).write_parquet(out/'candidates.parquet')
pl.DataFrame({'peptide_id':['run/100/2/50/10/1/ACDEFG'],
    'seq':['ACDEFG'], 'local_confidence':[[float('nan')]*6], 'search_type':['db'],
    'alc':[float('nan')], 'alc_rank':[1.], 'n_reads':[1]}).write_parquet(out/'peptides.parquet')
(out/'model.json').write_text(json.dumps(dict(conf_x=[0.,1.],conf_e=[.1,.01],
    substitution=mass_prior().tolist(), e_db=.1)))
(out/'samples.csv').write_text('sample,candidates,peptides,producer,model\n'+
    f'technical,{out}/candidates.parquet,{out}/peptides.parquet,technical fixture,{out}/model.json\n')
