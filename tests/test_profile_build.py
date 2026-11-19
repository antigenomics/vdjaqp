"""Resource accounting includes descendants once and excludes unrelated jobs."""
import importlib.util
import pytest
from pathlib import Path


def test_process_tree_and_owned_group_memory():
    path = Path(__file__).parents[1]/'scripts/profile_build.py'
    spec = importlib.util.spec_from_file_location('profile_build', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rows = [(10,1,10,100,3.), (20,10,10,200,80.), (30,20,30,300,20.),
            (40,1,10,400,0.), (50,1,50,999999,90.)]
    selected = module.process_scope(rows,10)
    assert {r[0] for r in selected} == {10,20,30,40}
    assert sum(r[3] for r in selected)*1024 == 1_024_000


def load_monitor():
    path = Path(__file__).parents[1]/'scripts/profile_build.py'
    spec = importlib.util.spec_from_file_location('profile_build_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_du_accepts_only_vanished_descendants_with_numeric_summary(tmp_path, monkeypatch):
    import subprocess
    import pytest
    module = load_monitor()
    root = tmp_path.resolve()
    seen = []
    def run(*args, **kwargs):
        seen.append(kwargs['env']['LC_ALL'])
        return subprocess.CompletedProcess(args[0], 1, f'17\t{root}\n', diagnostic)
    monkeypatch.setattr(module.subprocess, 'run', run)
    for diagnostic in [
        f"du: cannot access '{root}/reads.fastq': No such file or directory\n",
        f"du: cannot read directory '{root}/mmseqs_tmp/search': No such file or directory\n",
    ]:
        assert module.du_allocated(root) == (17*1024, True)
    assert seen == ['C', 'C']
    for status, stdout, stderr in [
        (1, f'17\t{root}\n', f'du: {root}/reads: Permission denied\n'),
        (1, f'17\t{root}\n', f"du: cannot read directory '{root}/reads': Permission denied\n"),
        (1, f'17\t{root}\n', f"du: cannot read directory '{root.parent}/outside': No such file or directory\n"),
        (1, 'malformed\n', f'du: {root}/reads: No such file or directory\n'),
        (1, f'17\t{root}\n', f'du: {root}: No such file or directory\n'),
        (1, f'17\t{root}\n', f'du: {root}/reads: No such file or directory\ndu: bad: Permission denied\n'),
        (2, f'17\t{root}\n', f'du: {root}/reads: No such file or directory\n'),
    ]:
        monkeypatch.setattr(module.subprocess, 'run', lambda *a, **k: subprocess.CompletedProcess(a[0], status, stdout, stderr))
        with pytest.raises(RuntimeError, match='du measurement failed'):
            module.du_allocated(root)


def test_monitor_exception_stops_owned_group_and_preserves_report(tmp_path, monkeypatch):
    import json
    import sys
    import pytest
    import vdjaqp
    module = load_monitor()
    output = tmp_path/'resource.json'
    workload = tmp_path/'workload'
    monkeypatch.setattr(sys, 'argv', ['profile_build', '--output', str(output),
        '--project-root', str(tmp_path), '--workload-root', str(workload), '--', 'test-command'])
    monkeypatch.setattr(vdjaqp, 'resolve_cpu_budget', lambda budget: 1)
    monkeypatch.setattr(module, 'du_allocated', lambda path: (100, False))
    def snapshot(pid):
        if pid == module.os.getpid(): return [], 0
        raise RuntimeError('exact mock monitoring failure')
    monkeypatch.setattr(module, 'snapshot', snapshot)
    class Child:
        pid = 987654
        returncode = None
        def poll(self): return self.returncode
    child = Child()
    monkeypatch.setattr(module.subprocess, 'Popen', lambda *a, **k: child)
    terminated = []
    def terminate(proc):
        terminated.append(proc.pid)
        proc.returncode = -15
    monkeypatch.setattr(module, 'terminate_owned', terminate)
    with pytest.raises(SystemExit) as exit_info:
        module.main()
    assert exit_info.value.code != 0
    report = json.loads(output.read_text())
    assert report['failure'] == 'monitor_error'
    assert report['monitor_error']['type'] == 'RuntimeError'
    assert report['monitor_error']['message'] == 'exact mock monitoring failure'
    assert 'RuntimeError: exact mock monitoring failure' in report['monitor_error']['traceback']
    assert report['exit_code'] == -15
    assert terminated == [child.pid]
    assert report['transient_partial_disk_scans'] == 0


@pytest.mark.parametrize("external", [False, True])
def test_workload_counts_roots_once_without_siblings(tmp_path, monkeypatch, external):
    import json
    import sys
    import pytest
    import vdjaqp
    module = load_monitor()
    project = tmp_path/'project'; project.mkdir()
    workload = (tmp_path/'external-work' if external else project)/'output'
    workload.mkdir(parents=True)
    sibling = tmp_path/'unrelated'; sibling.mkdir()
    output = workload.parent/'resource.json'
    measured = []
    def allocated(path):
        path = Path(path).resolve(); measured.append(path)
        assert path in (project, workload)
        return ((100 if external else 300) if path == project else 200), False
    monkeypatch.setattr(module, 'du_allocated', allocated)
    monkeypatch.setattr(module, 'snapshot', lambda pid: ([], 0))
    monkeypatch.setattr(vdjaqp, 'resolve_cpu_budget', lambda budget: 1)
    launched = []
    class Child:
        pid = 987654
        def poll(self): return 0
    def launch(*args, **kwargs):
        launched.append(kwargs['cwd']); return Child()
    monkeypatch.setattr(module.subprocess, 'Popen', launch)
    arguments = ['profile_build', '--output', str(output), '--project-root', str(project),
                 '--workload-root', str(workload), '--disk-limit', '400', '--', 'test-command']
    monkeypatch.setattr(sys, 'argv', arguments)
    module.main()
    report = json.loads(output.read_text())
    assert launched == [project]
    assert report['failure'] is None
    assert report['final_project_allocated_bytes'] == 300
    assert sibling not in measured and tmp_path not in measured
    output.unlink(); output.with_suffix('.log').unlink(); launched.clear()
    arguments[arguments.index('400')] = '250'
    with pytest.raises(SystemExit): module.main()
    assert launched == []


@pytest.mark.parametrize('timeout', ['nan', 'inf', '-inf'])
def test_nonfinite_timeout_rejected_before_monitoring(tmp_path, monkeypatch, timeout):
    import sys
    import vdjaqp
    module = load_monitor()
    output = tmp_path/'report.json'
    monkeypatch.setattr(sys, 'argv', ['profile_build', '--output', str(output),
        '--project-root', str(tmp_path), '--workload-root', str(tmp_path/'work'),
        '--timeout='+timeout, '--', 'unused-command'])
    monkeypatch.setattr(vdjaqp, 'resolve_cpu_budget', lambda budget: 1)
    monkeypatch.setattr(module, 'snapshot', lambda pid: pytest.fail('monitoring started'))
    with pytest.raises(SystemExit) as error:
        module.main()
    assert error.value.code == 2
    assert not output.exists()


def test_termination_stops_group_after_leader_exits(tmp_path):
    import os
    import signal
    import subprocess
    import sys
    import time
    module = load_monitor()
    marker = tmp_path/'child.pid'
    child = ('import os,signal,time;from pathlib import Path;'
        'signal.signal(signal.SIGTERM,signal.SIG_IGN);'
        f'Path({str(marker)!r}).write_text(str(os.getpid()));time.sleep(30)')
    parent = ('import subprocess,sys,time;'
        f'subprocess.Popen([sys.executable,"-c",{child!r}]);time.sleep(30)')
    proc = subprocess.Popen([sys.executable, '-c', parent], start_new_session=True)
    try:
        deadline = time.monotonic()+5
        while not marker.exists() and time.monotonic()<deadline:
            time.sleep(.02)
        assert marker.exists(), 'synthetic descendant did not start'
        pid = int(marker.read_text())
        module.terminate_owned(proc)
        # A reparented killed child can remain a zombie until the system reaps it.
        while True:
            state = subprocess.run(['ps','-p',str(pid),'-o','stat='],
                text=True,capture_output=True,timeout=2).stdout.strip()
            if not state or state.startswith('Z'):
                break
            if time.monotonic()>=deadline:
                pytest.fail('descendant survived owned-group termination')
            time.sleep(.02)
        assert proc.poll() is not None
    finally:
        try:
            os.killpg(proc.pid,signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait(timeout=5)
