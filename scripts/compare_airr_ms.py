"""Compare final AIRR-conditioned MS tables with recorded experimental targets.

No spectra, historical scoring or parameter fitting is run. Residuals retain each
target's units and denominators; recovery never enters the comparison objective.
"""
import argparse
import hashlib
import json
from pathlib import Path

import polars as pl


def target_residuals(targets, inventory, coverage, *, case, enzyme):
    """Compare panel-variable compatible species and matching residue unions only."""
    real = targets.filter(
        pl.col('condition').is_in(['solution_mixture','gel_mixture']) &
        (pl.col('reference')=='all') & (pl.col('enzyme')==enzyme) &
        pl.col('quantity').is_in(['compatible_called_species','residue_union_coverage']))
    expected={'FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4','junction','variable'}
    keys=['condition','quantity','region']
    if real.unique(keys).height!=real.height or coverage['region'].n_unique()!=coverage.height:
        raise ValueError('duplicate target or simulated region keys')
    if set(coverage['region'])!=expected:
        raise ValueError('expected all nine simulated panel regions')
    for condition in ('solution_mixture','gel_mixture'):
        group=real.filter(pl.col('condition')==condition)
        if (group.height!=10 or group.filter(pl.col('quantity')=='compatible_called_species').height!=1
            or set(group.filter(pl.col('quantity')=='residue_union_coverage')['region'])!=expected):
            raise ValueError('expected nine coverage targets and one species target per condition/enzyme')
    species = inventory.filter(pl.col('quantity')=='compatible_species')['count'].item()
    regions = {r['region']:r for r in coverage.to_dicts()}
    rows=[]
    for r in real.to_dicts():
        if r['quantity']=='compatible_called_species':
            value,numerator,denominator=species,species,None
        else:
            observed=regions[r['region']]
            if observed['region_length']!=r['denominator'] or observed['physical_reference_count']!=12:
                raise ValueError('simulation and recorded region populations/denominators differ')
            value,numerator,denominator=observed['coverage'],observed['covered_residues'],observed['region_length']
        if value is None:
            raise ValueError('a defined experimental target has undefined simulated coverage')
        rows.append(dict(case=case,condition=r['condition'],enzyme=enzyme,region=r['region'],
            quantity=r['quantity'],units=r['units'],real_value=r['value'],simulated_value=float(value),
            real_numerator=r['numerator'],real_denominator=r['denominator'],
            simulated_numerator=numerator,simulated_denominator=denominator,
            signed_residual=float(value)-r['value'],absolute_residual=abs(float(value)-r['value']),
            physical_reference_count=12,
            measurement_group='shared_chymotrypsin' if enzyme=='chymotrypsin' else r['condition']+'_'+enzyme,
            dependence='shared acquisition' if enzyme=='chymotrypsin' else 'pooled union; overlaps component acquisitions' if enzyme=='pooled' else 'separate acquisition',
            comparison='supplied heavy-variable panel; native search versus retained experimental exports'))
    return pl.DataFrame(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--targets',type=Path,required=True)
    parser.add_argument('--case',type=Path,action='append',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('output must be a fresh absent directory')
    targets=pl.read_csv(args.targets,separator='\t')
    tables=[];sources={}
    def digest(path):
        with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
    for case in args.case:
        report=case/'airr_ms_report.json'
        metadata=json.loads(report.read_text())
        if metadata['panel_rows']!=12:raise ValueError('expected the declared twelve-reference comparison')
        sources[str(report)]=digest(report)
        for enzyme,label in [('trypsin','Trypsin'),('chymotrypsin','Chymotrypsin'),('pooled','pooled')]:
            inv=case/(label+'_panel_species.parquet')
            cov=case/(label+'_panel_coverage_summary.parquet')
            sources[str(inv)]=digest(inv);sources[str(cov)]=digest(cov)
            tables.append(target_residuals(targets,pl.read_parquet(inv),pl.read_parquet(cov),case=case.name,enzyme=enzyme))
    result=pl.concat(tables,how='vertical_relaxed')
    args.output.mkdir(parents=True)
    result.write_csv(args.output/'target_residuals.tsv',separator='\t')
    result.write_parquet(args.output/'target_residuals.parquet',compression='zstd')
    (args.output/'comparison.json').write_text(json.dumps(dict(
        targets_sha256=digest(args.targets),source_sha256=sources,rows=result.height,
        selection='no weighted objective or automatic selection; yield/coverage residuals only; recovery excluded',
        scope='conditional heavy-variable comparison; total native species include constant stubs and are not matched experimental total yield',
        dependence='shared chymotrypsin counted once in any later joint objective; pooled unions are correlated'),indent=2)+'\n')


if __name__=='__main__':main()
