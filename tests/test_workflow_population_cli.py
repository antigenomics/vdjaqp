"""Invalid declared populations must fail before input loading or output creation."""
import subprocess
import sys
from pathlib import Path


def test_invalid_population_fails_before_loading(tmp_path):
    script=Path(__file__).resolve().parents[1]/'scripts/workflow_real.py'
    output=tmp_path/'output'
    result=subprocess.run([sys.executable,str(script),'--input-root',str(tmp_path/'missing'),'--output',str(output),'--expected-candidates','0'],capture_output=True,text=True)
    assert result.returncode==2
    assert 'expected candidate count must be positive' in result.stderr
    assert not output.exists()


def test_unknown_context_fails_before_loading(tmp_path):
    script=Path(__file__).resolve().parents[1]/'scripts/workflow_real.py'
    output=tmp_path/'output'
    result=subprocess.run([sys.executable,str(script),'--input-root',str(tmp_path/'missing'),'--output',str(output),'--case','unrecognised'],capture_output=True,text=True)
    assert result.returncode==2
    assert 'invalid choice' in result.stderr
    assert not output.exists()
