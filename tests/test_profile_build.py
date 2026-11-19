"""Resource accounting includes descendants once and excludes unrelated jobs."""
import importlib.util
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
