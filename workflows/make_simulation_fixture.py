"""Invent a tiny annotated protein substrate for native MS technical acceptance.

These are arbitrary protein strings with declared test coordinates, not a
biological repertoire, empirical observations or a simulation realism result.
"""
from pathlib import Path
import sys
import polars as pl

out=Path(sys.argv[1]).resolve()
out.mkdir(parents=True,exist_ok=False)
base='EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYY'
junction='CARGDYAMDYW'
suffix='GQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKS'
labels=[0]*20+[1]*10+[2]*15+[3]*10+[4]*(len(base)-55)+[6]*len(junction)+[7]*10+[8]*(len(suffix)-10)
rows=[]
for i in range(13):
    sequence=base+junction+suffix
    # Each candidate differs at one declared variable position.
    sequence=sequence[:i+2]+'A'+sequence[i+3:]
    rows.append(dict(clone_id=f'test{i:02}',seq=sequence,cdr3_seq=junction,
        cdr3_start=len(base),cdr3_end=len(base)+len(junction),region_ids=labels,
        v_name='IGHV1-2',j_name='IGHJ4',c_name='IGHG1',background_member=i>=12))
pl.DataFrame(rows).write_parquet(out/'candidates.parquet')
pl.DataFrame({'clone_id':[r['clone_id'] for r in rows[:12]]}).write_csv(out/'panel.csv')
