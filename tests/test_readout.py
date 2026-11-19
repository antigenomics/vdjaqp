"""Native spectrum identity must survive PSM adaptation and aggregation."""
import polars as pl
import pytest
from vdjaqp.aggregate import spectrum_key


def test_native_spectrum_grouping_and_missing_id():
    frame=pl.DataFrame({'peptide_id':['a','b'],'spectrum_id':['scan=5','scan=5']})
    assert spectrum_key(frame)['spectrum_id'].to_list()==['scan=5','scan=5']
    with pytest.raises(ValueError):
        spectrum_key(frame.with_columns(spectrum_id=pl.lit(None)))


def test_idxml_retains_native_metadata_and_missing_confidence(tmp_path):
    from vdjaqp.readout import read_idxml
    path=tmp_path/'psms.idXML'
    path.write_text('''<IdXML><IdentificationRun><PeptideIdentification RT="10" MZ="500" score_type="hyperscore"><UserParam name="spectrum_reference" value="controller/scan=9"/><PeptideHit sequence="ACDEFG" charge="2" score="33"/></PeptideIdentification></IdentificationRun></IdXML>''')
    table=read_idxml(path,run_id='synthetic')
    assert table['native_spectrum_reference'][0]=='controller/scan=9'
    assert table['hyperscore'][0]==33
    assert table['local_confidence'].list.len()[0]==6
    assert table['alc'].is_nan().all()
    assert spectrum_key(table)['spectrum_id'][0]==table['spectrum_id'][0]


@pytest.mark.parametrize("constant_supported", [False, True])
def test_observed_read_adapter_preserves_reverse_sequence_and_missing_count(tmp_path, constant_supported):
    import runpy
    from pathlib import Path
    from arda.refbuild.translate import reverse_complement
    export=runpy.run_path(str(Path(__file__).parents[1]/'scripts/smoke_repertoire.py'))['export_candidates']
    nt='GCTTGTGCTTGGTGGGGTGGTTCTTCT'
    junction=nt[3:15]
    rows=pl.DataFrame({'sequence_id':['read0'],'sequence':[reverse_complement(nt)],
        'rev_comp':['T'],'locus':['IGH'],'junction':[junction], 'junction_aa':['CAWW'],
        'v_call':['IGHV1-2'],'j_call':['IGHJ4'],'c_call':['IGHG1'],
        'fwr3_start':['1'],'fwr3_end':['6'],'cdr3_start':['7'],'cdr3_end':['12'],
        'fwr4_start':['13'],'fwr4_end':['21'],'mmseqs2_score':['100']})
    clones=pl.DataFrame({'locus':['IGH','IGH'],'junction':[junction,junction],
        'junction_aa':['CAWW','CAWW'],'v_call':['IGHV1-2','IGHV3-23'],'j_call':['IGHJ4']*2,
        'c_call':['IGHG1']*2,'duplicate_count':['3','1']})
    if constant_supported:
        rows=rows.with_columns(mmseqs2_tend=pl.lit('27'),mmseqs2_t_vjend=pl.lit('21'),mmseqs2_qend=pl.lit('27'))
    rows.write_csv(tmp_path/'sample.airr.tsv',separator='\t')
    clones.write_csv(tmp_path/'sample.clones.tsv',separator='\t')
    pl.DataFrame({'sequence_id':['read0'],'clone_row':[0]}).write_csv(tmp_path/'sample.read_ids.tsv',separator='\t')
    report=export(tmp_path/'sample',tmp_path/'observed.parquet')
    result=pl.read_parquet(tmp_path/'observed.parquet')
    assert result['seq'].to_list()==['ACAWWGGSS' if constant_supported else 'ACAWWGG','CAWW']
    assert result['coverage_kind'].to_list()==['observed_variable_fragment','junction_only']
    assert result['observed_read_id'][1] is None
    assert result['region_ids'][0].to_list()==[4,4,6,6,7,7,7]+([8,8] if constant_supported else [])
    assert result['junction_core_status'][0]=='conservative_whole_junction_unresolved'
    assert result['cdr3_start'][0]==1 and result['cdr3_end'][0]==5
    assert report['retained_observed_candidates']==2 and report['missing_observed_sequence']==0


def test_repertoire_cli_refuses_existing_stage(tmp_path):
    import subprocess
    import sys
    from pathlib import Path
    stage=tmp_path/'arda-amplicon';stage.mkdir()
    marker=stage/'keep.txt';marker.write_text('original')
    result=subprocess.run([sys.executable,str(Path(__file__).parents[1]/'scripts/smoke_repertoire.py'),'amplicon','--root',str(tmp_path)],capture_output=True,text=True)
    assert marker.exists() and marker.read_text()=='original'
    assert result.returncode!=0
    assert 'exists' in result.stderr


def test_native_junction_fallback_retains_partial_rows_and_counts_invalid(tmp_path):
    import runpy
    from pathlib import Path
    export=runpy.run_path(str(Path(__file__).parents[1]/'scripts/smoke_repertoire.py'))['export_candidates']
    pl.DataFrame(schema={'sequence_id':pl.String,'sequence':pl.String,'rev_comp':pl.String,
        'locus':pl.String,'junction':pl.String,'junction_aa':pl.String,'v_call':pl.String,'j_call':pl.String}).write_csv(tmp_path/'sample.airr.tsv',separator='\t')
    native=pl.DataFrame({'locus':['IGH']*3,'junction':['TGTGCTTGG','TGTGCTTGG','TGTAAATGG'],
        'junction_aa':['CAW','CCW','CAW'],'v_call':['IGHV1-2']*3,'j_call':['IGHJ4']*3,
        'c_call':['IGHG1']*3,'duplicate_count':['4','2','1']})
    pl.concat([native.head(1).with_columns(locus=pl.lit('IGK')),native]).write_csv(tmp_path/'sample.clones.tsv',separator='\t')
    pl.DataFrame(schema={'sequence_id':pl.String,'clone_row':pl.Int64}).write_csv(tmp_path/'sample.read_ids.tsv',separator='\t')
    report=export(tmp_path/'sample',tmp_path/'observed.parquet')
    result=pl.read_parquet(tmp_path/'observed.parquet')
    assert result.height==1
    assert result['seq'][0]=='CAW' and result['nt_seq'][0]=='TGTGCTTGG'
    assert result['coverage_kind'][0]=='junction_only'
    assert result['sequence_origin'][0]=='native_corrected_junction'
    assert result['native_clone_row'][0]==1
    assert result['region_ids'][0].to_list()==[4,6,7]
    assert result['junction_core_status'][0]=='conservative_whole_junction_unresolved'
    assert result['cdr3_start'][0]==0 and result['cdr3_end'][0]==3
    assert report['junction_only_candidates']==1 and report['invalid_native_junctions']==2
    assert report['missing_observed_sequence']==2


def test_fastp_overlap_preserves_unmerged_mates(tmp_path):
    import gzip,runpy,random
    from pathlib import Path
    repo=Path(__file__).parents[1]
    executable=repo/'.tool-envs/preprocess/bin/fastp'
    if not executable.exists():pytest.skip('isolated fastp runtime not installed')
    merge=runpy.run_path(str(repo/'scripts/smoke_preprocess.py'))['preprocess_pairs']
    # Exact software-string fixture, not a biological/read-error simulator.
    rng=random.Random(13);dna=lambda n:''.join(rng.choices('ACGT',k=n))
    template=dna(400);rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
    pairs=[(template[:300],rc(template[-300:])),(dna(300),dna(300))]
    for mate in [0,1]:
        with gzip.open(tmp_path/f'r{mate}.fq.gz','wt') as f:
            for i,pair in enumerate(pairs):f.write(f'@pair{i}/{mate+1}\n{pair[mate]}\n+\n'+('I'*300)+'\n')
    report=merge(tmp_path/'r0.fq.gz',tmp_path/'r1.fq.gz',tmp_path/'merged.fq.gz',threads=2,fastp=executable)
    assert report['input_pairs']==2 and report['merged_pairs']==1
    assert report['unmerged_pairs']==1 and report['retained_output_records']==3


def test_empty_native_search_export_remains_typed(tmp_path):
    from vdjaqp.readout import read_idxml
    path = tmp_path / 'empty.idXML'
    path.write_text('<IdXML><IdentificationRun/></IdXML>')
    frame = read_idxml(path, run_id='empty_native_search')
    assert frame.height == 0
    assert frame.schema['seq'] == pl.String
    assert frame.schema['local_confidence'] == pl.List(pl.Float32)
