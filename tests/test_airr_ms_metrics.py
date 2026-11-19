"""Supplied AIRR boundary and species checks; no native MS or cohort run."""
import importlib.util
from pathlib import Path
import polars as pl
import pytest
from vdjaqp.airr_ms_metrics import variable_sequences,supplied_reference_regions,eligible_species,species_metrics


def reference():
    # Explicit supplied region labels; junction C/W anchors belong to frameworks.
    return pl.DataFrame({'clone_id':['ref'],'seq':['AAAAACDEFGWTTSS'],
        'region_ids':[[0,1,2,3,4,5,5,6,6,5,5,7,7,8,8]],
        'cdr3_seq':['CDEFGW'],'cdr3_start':[5],'cdr3_end':[11]})


def test_supplied_boundaries_and_species_union():
    source=reference()
    assert variable_sequences(source)['variable_aa'].item()=='AAAAACDEFGWTT'
    regions=supplied_reference_regions(source)
    rows={r['region']:r for r in regions.to_dicts()}
    assert (rows['FR3']['start'],rows['FR3']['end'])==(4,6)
    assert (rows['CDR3']['start'],rows['CDR3']['end'])==(6,10)
    assert (rows['FR4']['start'],rows['FR4']['end'])==(10,13)
    calls=pl.DataFrame({'seq':['ACDEFI','ACDEFL','AAA','ACDEXI']})
    assert eligible_species(calls)['seq'].to_list()==['ACDEFL']
    queries=pl.DataFrame({'seq':['ACDEFL'],'peptide_id':['q']})
    hits=pl.DataFrame({'chain_id':['ref','ref'],'peptide_id':['q','q'],'start':[5,5],'end':[11,11]})
    inventory,coverage,total=species_metrics(calls,queries,hits,regions)
    counts=dict(zip(inventory['quantity'],inventory['count']))
    assert counts==dict(psm_assignment_rows=4,noncanonical_rows=1,short_rows=1,
                        eligible_species=1,compatible_species=1,unmapped_species=0)
    assert coverage.filter(pl.col('region')=='variable')['covered_residues'].item()==6
    assert total.filter(pl.col('region')=='variable')['region_length'].item()==13
    bad=source.with_columns(pl.Series('region_ids',[[8]+[0]*14]))
    with pytest.raises(ValueError,match='suffix'):variable_sequences(bad)
    bad=source.with_columns(pl.lit('CDEFGF').alias('cdr3_seq'))
    with pytest.raises(ValueError,match='junction'):supplied_reference_regions(bad)
    missing=source.with_columns(pl.Series('region_ids',[[0,1,2,3,4,5,5,6,6,5,5,0,0,8,8]]))
    with pytest.raises(ValueError,match='noncontiguous'):supplied_reference_regions(missing)


def test_driver_fresh_explicit_preflight_without_native(tmp_path):
    spec=importlib.util.spec_from_file_location('simulate_airr_ms',Path(__file__).parents[1]/'scripts/simulate_airr_ms.py')
    driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
    repertoire=tmp_path/'source.parquet';repertoire.touch()
    targets=tmp_path/'targets.csv';targets.touch()
    flags=['--repertoire',str(repertoire),'--target-ids',str(targets),'--source-format','supplied',
           '--output',str(tmp_path/'fresh'),'--physical-scope','panel','--ms2-per-bin','5']
    args=driver.parse_args(flags)
    assert args.physical_scope=='panel' and args.ms2_per_bin==5
    assert args.native_protein_intensity==10000. and args.stage_timeout==900.
    executable=tmp_path/'MSSimulator-fixed';executable.write_text('#!/bin/sh\nexit 0\n')
    with pytest.raises(SystemExit):driver.parse_args(flags+['--mssimulator',str(executable)])
    executable.chmod(0o755)
    assert driver.parse_args(flags+['--mssimulator',str(executable)]).mssimulator==executable.resolve()
    for timeout in ('0','nan','inf'):
        with pytest.raises(SystemExit):driver.parse_args(flags+['--stage-timeout',timeout])
    for intensity in ('0','-1','nan','inf'):
        with pytest.raises(SystemExit):driver.parse_args(flags+['--native-protein-intensity',intensity])
    with pytest.raises(SystemExit):driver.parse_args(flags+['--ms2-per-bin','7'])
    (tmp_path/'fresh').mkdir()
    with pytest.raises(SystemExit):driver.parse_args(flags)


def test_driver_orchestration_keeps_search_source_and_global_aliases(tmp_path,monkeypatch):
    import json
    import vdjaqp.pipeline
    import vdjaqp.readout
    spec=importlib.util.spec_from_file_location('airr_driver',Path(__file__).parents[1]/'scripts/simulate_airr_ms.py')
    driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
    source=pl.concat([reference().with_columns(pl.lit(f'id{i}').alias('clone_id')) for i in range(12)])
    source=source.with_columns(pl.lit('IGHV1-2*01').alias('v_name'),pl.lit('IGHJ1*01').alias('j_name'),pl.lit('IGHG1*01').alias('c_name'))
    repertoire=tmp_path/'input.parquet';source.write_parquet(repertoire)
    targets=tmp_path/'targets.csv';source.select('clone_id').write_csv(targets)
    root=tmp_path/'fresh'
    executable=tmp_path/'MSSimulator-fixed';executable.write_text('#!/bin/sh\nexit 0\n');executable.chmod(0o755)
    calls=[]
    def fake_native(command,**kwargs):
        out=Path(command[4]);out.mkdir()
        (out/'ms2_truth.json').write_text('[]')
        (out/'spectra.mzML').write_text('synthetic test spectra placeholder')
        calls.append(command)
    monkeypatch.setattr(driver.subprocess,'run',fake_native)
    monkeypatch.setattr(vdjaqp.readout,'read_idxml',lambda *a,**kw:pl.DataFrame(schema=vdjaqp.readout._SCHEMA))
    class FakeEvidence:
        @classmethod
        def from_frames(cls,candidates,peptides,**kwargs):
            assert candidates['v_name'].unique().to_list()==['IGHV1-2']
            assert candidates['seq'].to_list()==source['seq'].to_list()
            obj=cls();obj.ids=candidates['clone_id'];return obj
        def score(self):
            return pl.DataFrame({'clone_id':self.ids,'score':[0.]*12,'confirmed_nats':[0.]*12})
    monkeypatch.setattr(vdjaqp.pipeline,'Evidence',FakeEvidence)
    driver.main(['--repertoire',str(repertoire),'--target-ids',str(targets),
        '--source-format','synthetic unit fixture','--output',str(root),'--physical-scope','all',
        '--ms2-per-bin','5','--cpus','1','--stage-timeout','901','--mssimulator',str(executable)])
    assert len(calls)==2 and all('1200' in command for command in calls)
    assert all(command[command.index('--stage-timeout')+1]=='901.0' for command in calls)
    assert all(command[command.index('--mssimulator')+1]==str(executable.resolve()) for command in calls)
    assert not list(root.rglob('spectra.mzML'))
    report=json.loads((root/'airr_ms_report.json').read_text())
    assert report['physical_rows']==report['candidate_rows']==12
    recovery=pl.read_parquet(root/'pooled_identity_metrics.parquet')
    assert recovery.filter(pl.col('quantity')=='recovered_classes_at20')['count'].item()==0
    assert report['readouts']['pooled']['panel_recovery']['variable_mass_classes']==1


def test_constant_stub_and_cross_boundary_hits_are_excluded_from_fit():
    from vdjaqp.airr_ms_metrics import supplied_raw_regions,variable_only_placements,supplied_gene_coverage
    source=reference().with_columns(pl.lit('IGHV1-2*01').alias('v_name'),
        pl.lit('IGHJ1*01').alias('j_name'),pl.lit('IGHG1*01').alias('c_name'))
    regions=supplied_raw_regions(source)
    raw={r['region']:r for r in regions.to_dicts()}
    assert (raw['Cstub']['start'],raw['Cstub']['end'])==(13,15)
    assert raw['Cstub']['gene_call']=='IGHG1'
    assert raw['rawpart5_run1']['end']==7 and raw['rawpart5_run2']['start']==9
    placements=pl.DataFrame({'chain_id':['ref','ref'],'peptide_id':['cross','within'],
        'start':[9,5],'end':[15,11]})
    metadata=pl.DataFrame({'chain_id':['ref'],'variable_length':[13]})
    fit_hits=variable_only_placements(placements,metadata)
    assert fit_hits['peptide_id'].to_list()==['within']
    queries=pl.DataFrame({'peptide_id':['cross','within'],'seq':['GWTTSI','ACDEFL']})
    calls=pl.DataFrame({'seq':['GWTTSI']})
    # I/L equivalence has already normalized the query lookup.
    queries=queries.with_columns(pl.col('seq').str.replace_all('I','L',literal=True))
    _,coverage,_=species_metrics(calls,queries,placements,regions)
    assert coverage.filter(pl.col('region')=='Cstub')['covered_residues'].item()==2
    gene=supplied_gene_coverage(coverage).filter(pl.col('region')=='Cstub').row(0,named=True)
    assert gene['gene_call']=='IGHG1' and gene['covered_residues']==gene['residue_denominator']==2
    inventory,panel_coverage,_=species_metrics(calls,queries,fit_hits,supplied_reference_regions(source))
    assert inventory.filter(pl.col('quantity')=='compatible_species')['count'].item()==0
    assert panel_coverage.filter(pl.col('region')=='variable')['covered_residues'].item()==0
