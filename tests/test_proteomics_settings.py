"""Reject invalid acquisition settings before native tools or output creation."""
from pathlib import Path
import subprocess
import sys
import pytest


def test_invalid_acquisition_rejected_before_output(tmp_path):
    script = Path(__file__).parents[1]/'scripts/smoke_proteomics.py'
    for flags in (['--gradient-seconds','nan'],['--ms2-per-bin','0']):
        out=tmp_path/'not-created'
        result=subprocess.run([sys.executable,str(script),'missing-proteins','missing-db',
                               str(out),*flags],capture_output=True,text=True,timeout=10)
        assert result.returncode != 0 and not out.exists()
        assert 'finite positive gradient and positive MS2 capacity' in result.stderr


def test_stale_native_manifest_rejected_before_output(tmp_path):
    script=Path(__file__).parents[1]/'scripts/smoke_proteomics.py'
    binary=tmp_path/'native-tool'
    binary.write_text('#!/bin/sh\nexit 99\n');binary.chmod(0o700)
    Path(str(binary)+'.manifest.json').write_text('{"executable_sha256":"stale"}')
    out=tmp_path/'not-created'
    result=subprocess.run([sys.executable,str(script),'missing-proteins','missing-db',str(out),
        '--mssimulator',str(binary)],capture_output=True,text=True,timeout=10)
    assert result.returncode!=0 and not out.exists()
    assert 'does not match its build manifest' in result.stderr


def test_executable_identity_records_actual_bytes_and_path_without_platform_claim(tmp_path):
    import hashlib
    import importlib.util
    script=Path(__file__).parents[1]/'scripts/smoke_proteomics.py'
    spec=importlib.util.spec_from_file_location('proteomics_identity',script)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for platform, signature in [('linux',b'\x7fELF'),('darwin',b'\xcf\xfa\xed\xfe')]:
        directory=tmp_path/platform
        directory.mkdir()
        for name in ('SimpleSearchEngine','TextExporter'):
            binary=directory/name
            content=signature+name.encode()
            binary.write_bytes(content)
            identity=module.executable_identity(binary)
            assert identity=={'path':str(binary.resolve()),'sha256':hashlib.sha256(content).hexdigest()}
            assert set(identity)=={'path','sha256'}


@pytest.mark.parametrize('flags', [
    ['--min-detectability','nan'], ['--min-detectability','inf'],
    ['--min-detectability','-0.01'], ['--min-detectability','1.01'],
    ['--dynamic-exclusion-seconds','nan'], ['--dynamic-exclusion-seconds','inf'],
    ['--dynamic-exclusion-seconds','-1']])
def test_invalid_native_assay_knobs_rejected_before_output(tmp_path, flags):
    script=Path(__file__).parents[1]/'scripts/smoke_proteomics.py'
    out=tmp_path/'not-created'
    result=subprocess.run([sys.executable,str(script),'missing','missing',str(out),*flags],
                          capture_output=True,text=True,timeout=10)
    assert result.returncode!=0 and not out.exists()
    assert 'must be finite' in result.stderr


@pytest.mark.parametrize('detectability, exclusion', [(.5,0),(.5,30),(.75,0),(.75,30)])
def test_native_assay_ini_values_and_default_identity(tmp_path, monkeypatch, detectability, exclusion):
    """Exercise the real CLI/INI-writing path with no native acquisition."""
    import importlib.util,json,xml.etree.ElementTree as ET
    script=Path(__file__).parents[1]/'scripts/smoke_proteomics.py'
    spec=importlib.util.spec_from_file_location('proteomics_settings',script)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    # Published native values already used by the existing baseline.
    baseline={
        'algorithm:MSSim:Digestion:enzyme':'Trypsin',
        'algorithm:MSSim:RT:total_gradient_time':'1200.0',
        'algorithm:MSSim:RT:scan_window:min':'0.0',
        'algorithm:MSSim:RT:scan_window:max':'1200.0',
        'algorithm:MSSim:RT:sampling_rate':'2.0',
        'algorithm:MSSim:RawSignal:enabled':'true',
        'algorithm:MSSim:RawTandemSignal:status':'precursor',
        'algorithm:MSSim:RawTandemSignal:tandem_mode':'0',
        'algorithm:MSSim:RawTandemSignal:Precursor:ms2_spectra_per_rt_bin':'5',
        'algorithm:MSSim:Detectability:dt_simulation_on':'true',
        'algorithm:MSSim:Detectability:min_detect':'0.5',
        'algorithm:MSSim:RawTandemSignal:Precursor:Exclusion:use_dynamic_exclusion':'false',
        'algorithm:MSSim:RawTandemSignal:Precursor:Exclusion:exclusion_time':'100.0',
        'algorithm:RandomNumberGenerators:biological':'reproducible',
        'algorithm:RandomNumberGenerators:technical':'reproducible'}
    def execute(command,**kwargs):
        if '-write_ini' in command:
            root=ET.Element('PARAMETERS');top=ET.SubElement(root,'NODE',name=Path(command[0]).name)
            node=ET.SubElement(top,'NODE',name='1')
            if Path(command[0]).name=='MSSimulator':
                for key,value in baseline.items():
                    parent=node
                    for part in key.split(':')[:-1]:
                        match=next((n for n in parent if n.get('name')==part),None)
                        parent=match if match is not None else ET.SubElement(parent,'NODE',name=part)
                    ET.SubElement(parent,'ITEM',name=key.split(':')[-1],value=value)
            ET.ElementTree(root).write(command[command.index('-write_ini')+1])
        elif '-out' in command:
            path=Path(command[command.index('-out')+1]);path.write_text('<root/>')
    monkeypatch.setattr(module.subprocess,'run',execute)
    monkeypatch.setattr(module,'executable_identity',lambda p:dict(path=str(p),sha256='fixture'))
    monkeypatch.setattr(module,'link_native_truth',lambda *a:([],set()))
    out=tmp_path/'draw'
    fake_tool=tmp_path/'MSSimulator';fake_tool.write_text('#!/bin/sh\nexit 0\n');fake_tool.chmod(0o700)
    argv=[str(script),'proteins','search',str(out),'--threads','1','--gradient-seconds','1200',
          '--mssimulator',str(tmp_path/'MSSimulator')]
    if (detectability,exclusion)!=(.5,0):argv+=['--min-detectability',str(detectability),'--dynamic-exclusion-seconds',str(exclusion)]
    monkeypatch.setattr(sys,'argv',argv);module.main()
    actual={}
    def flatten(node,parts=()):
        for child in node:
            path=(*parts,child.get('name',''))
            if child.tag=='ITEM':actual[':'.join(path[2:])]=child.get('value')
            flatten(child,path)
    flatten(ET.parse(out/'mssimulator.ini').getroot())
    expected=dict(baseline);expected['algorithm:MSSim:Detectability:min_detect']=str(detectability)
    if exclusion:
        expected['algorithm:MSSim:RawTandemSignal:Precursor:Exclusion:use_dynamic_exclusion']='true'
        expected['algorithm:MSSim:RawTandemSignal:Precursor:Exclusion:exclusion_time']=str(float(exclusion))
    assert actual==expected
    report=json.loads((out/'summary.json').read_text())
    assert report['min_detectability']==detectability and report['dynamic_exclusion_seconds']==exclusion
