"""Fixed-input abundance bookkeeping and 255-target selection; no native draws."""
import copy,json,math
import numpy as np
import polars as pl
import pytest
from test_adapter import example
from vdjaqp.airr_ms_inputs import prepare_airr_ms_inputs
from vdjaqp.abundance_inference import relative_amounts,target_vector,distance,REGIONS,BLOCKS


def targets():
    rows=[]
    for condition,enzyme in [*BLOCKS,('gel_mixture','chymotrypsin')]:
        rows.append(dict(condition=condition,enzyme=enzyme,reference='all',region='',quantity='compatible_called_species',
                         value=4.,numerator=4,denominator=None,units='species'))
        for i in range(1,13):
            for region in REGIONS:
                rows.append(dict(condition=condition,enzyme=enzyme,reference=f'ref{i}',region=region,
                    quantity='residue_union_coverage',value=.5,numerator=5,denominator=10,units='fraction'))
    return pl.DataFrame(rows)


def test_prior_summary_and_distance_exact_contract():
    assert np.array_equal(relative_amounts(np.zeros(11)),np.ones(12))
    assert relative_amounts(np.arange(11.)).sum()==pytest.approx(12.)
    t=target_vector(targets());assert np.asarray(t['coverage']).size+len(t['yield_counts'])==255
    assert distance(t,t)==0
    changed=copy.deepcopy(t);changed['coverage'][0][0][0]=1.
    assert distance(changed,t)==pytest.approx(.5*.5/(84*3))
    changed['denominators'][0][0][0]=11
    with pytest.raises(ValueError,match='denominators'):distance(changed,t)
    with pytest.raises(ValueError):target_vector(pl.concat([targets(),targets().head(1)]))
    with pytest.raises(ValueError):relative_amounts([float('nan')]*11)


def test_explicit_relative_amounts_and_unchanged_sequence_inputs(tmp_path):
    frame=pl.concat([example().with_columns(clone_id=pl.lit(f'c{i}')) for i in range(12)])
    amounts={f'c{i}':.5 if i<6 else 1.5 for i in range(12)}
    paths=prepare_airr_ms_inputs(frame,physical_clone_ids=list(amounts),equal_protein_amounts=False,
        protein_relative_amounts=amounts,source_format='synthetic fixed export',output=tmp_path/'weighted')
    assert paths['proteins_fasta'].read_text().count('intensity=5000')==6
    assert paths['proteins_fasta'].read_text().count('intensity=15000')==6
    manifest=json.loads(paths['manifest'].read_text());assert manifest['total_native_input_intensity']==120000
    assert pl.read_parquet(paths['candidates']).equals(frame)
    assert pl.read_csv(paths['id_mapping'],separator='\t')['protein_relative_amount'].sum()==12
    for change in [{**amounts,'extra':1},dict(amounts,c0=0),dict(amounts,c0=float('nan')),dict(amounts,c0=1)]:
        with pytest.raises(ValueError):
            prepare_airr_ms_inputs(frame,physical_clone_ids=list(amounts),equal_protein_amounts=False,
                protein_relative_amounts=change,source_format='synthetic fixed export',output=tmp_path/'invalid')
        assert not (tmp_path/'invalid').exists()


@pytest.mark.parametrize('assay', [{},dict(min_detectability=.75,dynamic_exclusion_seconds=30.)])
def test_joint_draw_has_three_fresh_acquisitions_and_one_shared_chymotrypsin(tmp_path,monkeypatch,assay):
    import vdjaqp.abundance_forward as joint
    import vdjaqp.pipeline
    from vdjaqp.readout import _SCHEMA
    from test_airr_ms_metrics import reference
    ref=reference().with_columns(v_name=pl.lit('IGHV1'),j_name=pl.lit('IGHJ1'),c_name=pl.lit('IGHG1'))
    candidates=pl.concat([ref.with_columns(clone_id=pl.lit(f'c{i}')) for i in range(33052)])
    repertoire=tmp_path/'input.parquet';candidates.write_parquet(repertoire)
    ids=tmp_path/'ids.csv';pl.DataFrame({'clone_id':[f'c{i}' for i in range(12)]}).write_csv(ids)
    native=[]
    def execute(command,**kwargs):
        native.append(command);out=__import__('pathlib').Path(command[4]);out.mkdir()
        (out/'ms2_truth.json').write_text('[]');(out/'spectra.mzML').write_bytes(b'fixture')
        (out/'psms.idXML').write_text('fixture');(out/'truth.featureXML').write_text('consumed fixture')
        (out/'summary.json').write_text(json.dumps({'fixture':True}))
    monkeypatch.setattr(joint.subprocess,'run',execute)
    monkeypatch.setattr(joint,'read_idxml',lambda *a,**k:pl.DataFrame(schema=_SCHEMA))
    monkeypatch.setattr(joint,'exact_peptide_placements',lambda *a,**k:pl.DataFrame(schema={
        'chain_id':pl.String,'peptide_id':pl.String,'start':pl.Int32,'end':pl.Int32}))
    scored=[]
    class Evidence:
        @classmethod
        def from_frames(cls,frame,peptides,**kwargs):
            assert frame['background_member'].sum()==33040
            assert frame.filter(~pl.col('background_member')).height==12
            scored.append(peptides.height);result=cls();result.ids=frame['clone_id'];return result
        def score(self):return pl.DataFrame({'clone_id':self.ids,'score':[0.]*len(self.ids)})
    monkeypatch.setattr(vdjaqp.pipeline,'Evidence',Evidence)
    out=tmp_path/'draw'
    summary=joint.forward_draw(repertoire,ids,'synthetic software fixture',np.ones(12),out,cpus=1,**assay)
    assert len(native)==3 and len(scored)==2
    assert sum(command[command.index('--enzyme')+1]=='Chymotrypsin' for command in native)==1
    assert summary['yield_counts']==[0,0,0]
    assert np.asarray(summary['coverage']).shape==(3,12,7)
    report=json.loads((out/'joint_draw.json').read_text())
    assert report['status']=='complete' and report['shared_chymotrypsin_acquisitions']==1
    assert all(v['status']=='deleted' for v in report['cleanup'].values())
    assert all(v['status']=='deleted' for v in report['consumed_native_features'].values())
    assert report['recovery_used_for_selection'] is False
    expected=dict(min_detectability=.5,dynamic_exclusion_seconds=0.)|assay
    assert report['assay_parameters']==expected
    for command in native:
        assert float(command[command.index('--min-detectability')+1])==expected['min_detectability']
        assert float(command[command.index('--dynamic-exclusion-seconds')+1])==expected['dynamic_exclusion_seconds']
    assert not list(out.rglob('spectra.mzML'))


@pytest.mark.parametrize('assay', [dict(min_detectability=x) for x in [-.01,1.01,float('nan'),float('inf')]]+
    [dict(dynamic_exclusion_seconds=x) for x in [-1.,float('nan'),float('inf')]])
def test_invalid_assay_settings_rejected_before_inputs_or_output(tmp_path,assay):
    from vdjaqp.abundance_forward import forward_draw
    out=tmp_path/'not-created'
    with pytest.raises(ValueError,match='must be finite'):
        forward_draw('missing','missing','software fixture',np.ones(12),out,cpus=1,**assay)
    assert not out.exists()


@pytest.mark.parametrize('flags, expected', [([],dict(min_detectability=.5,dynamic_exclusion_seconds=0.)),
    (['--min-detectability','.75','--dynamic-exclusion-seconds','30'],
     dict(min_detectability=.75,dynamic_exclusion_seconds=30.))])
def test_forward_cli_propagates_assay_settings(tmp_path,monkeypatch,flags,expected):
    import runpy,sys
    from pathlib import Path
    import vdjaqp.abundance_forward as joint
    observed=target_vector(targets());received=[]
    def forward(*args,**kwargs):
        received.append(kwargs);Path(args[4]).mkdir();return observed
    monkeypatch.setattr(joint,'forward_draw',forward)
    target=tmp_path/'targets.tsv';targets().write_csv(target,separator='\t')
    amounts=tmp_path/'amounts.json';amounts.write_text(json.dumps([1.]*12))
    output=tmp_path/'out';script=Path(__file__).parents[1]/'scripts/run_abundance_forward.py'
    monkeypatch.setattr(sys,'argv',[str(script),'--repertoire','missing','--target-ids','missing',
        '--source-format','software fixture','--output',str(output),'--targets',str(target),
        '--amounts',str(amounts),'--cpus','1',*flags])
    runpy.run_path(str(script),run_name='__main__')
    assert len(received)==1
    assert {k:received[0][k] for k in expected}==expected
    discrepancy=json.loads((output/'discrepancy.json').read_text())
    assert discrepancy['distance']==0 and discrepancy['coordinates']==255
    assert discrepancy['assay_parameters']==expected and discrepancy['recovery_used'] is False


def test_installed_published_engine_counts_prior_calibration_and_hard_guard(tmp_path):
    import os,subprocess
    from pathlib import Path
    repo=Path(__file__).parents[1];python=repo/'.tool-envs/inference/bin/python'
    if not python.exists():pytest.skip('qualified isolated pyABC environment required')
    code='''
import tempfile
from pathlib import Path
from vdjaqp.abundance_inference import run_abc,target_vector
rows=[]
for condition,enzyme in [('solution_mixture','trypsin'),('gel_mixture','trypsin'),('solution_mixture','chymotrypsin'),('gel_mixture','chymotrypsin')]:
 rows.append(dict(condition=condition,enzyme=enzyme,reference='all',region='',quantity='compatible_called_species',value=1,numerator=1,denominator=None,units='species'))
 for i in range(1,13):
  for region in ['FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4']:
   rows.append(dict(condition=condition,enzyme=enzyme,reference=f'ref{i}',region=region,quantity='residue_union_coverage',value=0,numerator=0,denominator=1,units='fraction'))
summary=target_vector(rows)  # Actual CLI uses list-of-dict CSV rows, not a DataFrame.
with tempfile.TemporaryDirectory() as temp:
 output=Path(temp)/'run';calls=[]
 report=run_abc(lambda amounts,index,seconds:(calls.append(index) or summary),summary,output,max_draws=100,wall_seconds=60)
 assert report['status']=='complete' and report['attempted_draws']==len(calls)==64
 assert report['completed_populations']==1 and report['epsilon_history']=={'0':.1}
 assert report['population_effective_sample_size']=={'0':32.}
 assert len(list(output.glob('attempt-*.json')))==64
with tempfile.TemporaryDirectory() as temp:
 calls=[];output=Path(temp)/'run'
 report=run_abc(lambda amounts,index,seconds:(calls.append(index) or summary),summary,output,max_draws=3,wall_seconds=60)
 assert report['status']=='budget_limited' and report['attempted_draws']==len(calls)==3
 assert report['completed_populations']==0 and report['epsilon_history']=={}
 assert len(list(output.glob('attempt-*.json')))==3
'''
    subprocess.run([str(python),'-c',code],cwd=repo,env=dict(os.environ,
        MPLCONFIGDIR=str(repo/'.tool-envs/inference/mpl'),XDG_CACHE_HOME=str(repo/'.tool-envs/inference/cache')),
        check=True,capture_output=True,text=True,timeout=60)
