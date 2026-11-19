"""Final-table residual comparison preserves units and population denominators."""
import importlib.util
from pathlib import Path
import polars as pl
import pytest

spec=importlib.util.spec_from_file_location('comparison',Path(__file__).parents[1]/'scripts/compare_airr_ms.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_residual_units_and_denominators():
    regions=['FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4','junction','variable']
    target=pl.concat([pl.DataFrame(dict(condition=[condition]*10,reference=['all']*10,
        enzyme=['chymotrypsin']*10,region=[None,*regions],
        quantity=['compatible_called_species',*['residue_union_coverage']*9],units=['species',*['fraction']*9],
        value=[5.,*[.5]*9],numerator=[5,*[10]*9],denominator=[None,*[20]*9]))
        for condition in ('solution_mixture','gel_mixture')])
    inventory=pl.DataFrame({'quantity':['compatible_species'],'count':[6]})
    coverage=pl.DataFrame(dict(region=regions,region_length=[20]*9,physical_reference_count=[12]*9,
        coverage=[.75]*9,covered_residues=[15]*9))
    result=module.target_residuals(target,inventory,coverage,case='panel5',enzyme='chymotrypsin')
    assert result['signed_residual'].to_list()==[1.,*[.25]*9]*2
    assert result['units'].to_list()==['species',*['fraction']*9]*2
    assert result['simulated_denominator'].to_list()==[None,*[20]*9]*2
    assert result['measurement_group'].to_list()==['shared_chymotrypsin']*20
    with pytest.raises(ValueError,match='denominators'):
        module.target_residuals(target,inventory,coverage.with_columns(region_length=pl.lit(21)),case='panel5',enzyme='chymotrypsin')
    for bad in (target.head(19),pl.concat([target,target.head(1)])):
        with pytest.raises(ValueError):
            module.target_residuals(bad,inventory,coverage,case='panel5',enzyme='chymotrypsin')
