"""Run one explicitly selected fresh joint draw under the baseline scientific runtime."""
import argparse,json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from vdjaqp import resolve_cpu_budget

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repertoire',type=Path,required=True);p.add_argument('--target-ids',type=Path,required=True)
    p.add_argument('--source-format',required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--amounts',type=Path,required=True,help='JSON list of twelve relative amounts, sum twelve')
    p.add_argument('--targets',type=Path,required=True)
    p.add_argument('--cpus',type=int,default=16);p.add_argument('--mssimulator',type=Path)
    p.add_argument('--min-detectability',type=float,default=.5)
    p.add_argument('--dynamic-exclusion-seconds',type=float,default=0.)
    a=p.parse_args();budget=resolve_cpu_budget(a.cpus)
    import polars as pl
    from vdjaqp.abundance_inference import target_vector,distance
    observed=target_vector(pl.read_csv(a.targets,separator='\t'))
    from vdjaqp.abundance_forward import forward_draw
    result=forward_draw(a.repertoire,a.target_ids,a.source_format,json.loads(a.amounts.read_text()),a.output,
                 cpus=budget,mssimulator=a.mssimulator,expected_denominators=observed['denominators'],
                 min_detectability=a.min_detectability,dynamic_exclusion_seconds=a.dynamic_exclusion_seconds)
    (a.output/'discrepancy.json').write_text(json.dumps(dict(distance=distance(result,observed),
        coordinates=255,recovery_used=False,
        assay_parameters=dict(min_detectability=a.min_detectability,dynamic_exclusion_seconds=a.dynamic_exclusion_seconds)),indent=2)+'\n')
