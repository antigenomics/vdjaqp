"""The CPU budget is a ceiling across the whole process tree."""

import os
import subprocess
import sys

import pytest

from vdjaqp.resources import DEFAULT_CPU_BUDGET, resolve_cpu_budget


@pytest.fixture(autouse=True)
def _restore_budget():
    """The budget is process-global; put it back so test order cannot matter."""
    yield
    resolve_cpu_budget()


def test_default_uses_all_available_cores():
    assert resolve_cpu_budget() == resolve_cpu_budget("all")


def test_serial_means_one_thread():
    assert resolve_cpu_budget("serial") == 1
    assert os.environ["POLARS_MAX_THREADS"] == "1"


@pytest.mark.parametrize("bad", [0, -4, "many", "0"])
def test_invalid_budget_is_rejected_not_clamped(bad):
    with pytest.raises(ValueError):
        resolve_cpu_budget(bad)


def test_request_beyond_allocation_is_clamped():
    assert resolve_cpu_budget(100_000) == resolve_cpu_budget("all")


def test_limits_reach_the_native_pools():
    """Environment values alone are not evidence; polars must actually resize."""
    code = (
        "import vdjaqp; vdjaqp.resolve_cpu_budget(2); "
        "import polars as pl; print(pl.thread_pool_size())"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert out.stdout.strip() == "2"


def test_polars_io_pool_is_bounded_after_query():
    code = (
        "import vdjaqp; vdjaqp.resolve_cpu_budget(2); import os, polars as pl; "
        "pl.DataFrame({'v':[1,2]}).select(pl.col('v').sum()); "
        "print(pl.thread_pool_size(), os.environ['POLARS_ASYNC_THREAD_COUNT'])"
    )
    out = subprocess.run([sys.executable, '-c', code], capture_output=True,
                         text=True, check=True, timeout=30)
    assert out.stdout.strip() == '2 1'


@pytest.mark.skipif(not sys.platform.startswith('linux'), reason='Linux OpenBLAS pool regression')
def test_numpy_does_not_create_an_auxiliary_compute_team():
    code = (
        "import vdjaqp; vdjaqp.resolve_cpu_budget(4); import numpy; "
        "from pathlib import Path; print(len(list(Path('/proc/self/task').iterdir())))"
    )
    out = subprocess.run([sys.executable, '-c', code], capture_output=True,
                         text=True, check=True, timeout=30)
    assert int(out.stdout.strip()) <= 2


def test_allocation_uses_existing_cgroup_affinity_utility_and_slurm(monkeypatch):
    import vdjtools.cores
    monkeypatch.setattr(vdjtools.cores, 'available_cores', lambda: 12)
    monkeypatch.setenv('SLURM_CPUS_PER_TASK', '3')
    assert resolve_cpu_budget() == 3
    assert resolve_cpu_budget(100) == 3
    assert resolve_cpu_budget(2) == 2
    monkeypatch.setenv('SLURM_CPUS_PER_TASK', '24')
    assert resolve_cpu_budget() == 12


@pytest.mark.parametrize('bad', [True, False, 1.5, 2.0, float('inf')])
def test_noninteger_budget_is_rejected(bad):
    with pytest.raises(ValueError):
        resolve_cpu_budget(bad)


def test_per_call_threads_are_bounded_without_mutating_budget():
    from vdjaqp.resources import bounded_threads
    resolve_cpu_budget(2)
    selected = resolve_cpu_budget(2)
    assert bounded_threads() == bounded_threads(0) == selected
    assert bounded_threads(100_000) == selected
    assert bounded_threads(1) == 1
    assert os.environ['POLARS_MAX_THREADS'] == str(selected)


@pytest.mark.parametrize('bad', [-1, True, False, 1.5, '2'])
def test_bad_per_call_threads_rejected(bad):
    from vdjaqp.resources import bounded_threads
    with pytest.raises(ValueError):
        bounded_threads(bad)


@pytest.mark.parametrize('bad', [-1, True, 1.5])
def test_native_callers_validate_threads_before_processing(bad):
    import polars as pl
    from vdjaqp.match import PeptideMatcher
    from vdjaqp.pgen import junction_pgen, insertion_pgen
    from vdjaqp.significance import peptide_evalues
    matcher = PeptideMatcher(['c0'], ['ACDEFGHKLMNP'])
    calls = [lambda: matcher.search([], [], threads=bad),
             lambda: junction_pgen([], threads=bad),
             lambda: insertion_pgen(pl.DataFrame(), threads=bad),
             lambda: peptide_evalues(pl.DataFrame(), pl.DataFrame(), threads=bad)]
    for call in calls:
        with pytest.raises(ValueError, match='threads'):
            call()


def test_pgen_native_call_cannot_exceed_selected_budget(monkeypatch):
    import numpy as np
    from vdjtools.model import native, bundled
    from vdjaqp.pgen import junction_pgen
    seen = []
    monkeypatch.setattr(bundled, 'load_bundled', lambda *a, **kw: object())
    def batch(model, sequences, threads):
        seen.append(threads)
        return np.ones(len(sequences))
    monkeypatch.setattr(native, 'pgen_aa_batch', batch)
    resolve_cpu_budget(1)
    junction_pgen(['CARGGDYW'], threads=100_000)
    assert seen == [1]
    assert os.environ['POLARS_MAX_THREADS'] == '1'


def test_late_smaller_budget_cannot_misreport_initialized_polars_pool():
    code = '''
import os
import vdjaqp
if vdjaqp.resolve_cpu_budget('all') < 2:
    raise SystemExit(77)
vdjaqp.resolve_cpu_budget(2)
import polars as pl
assert pl.thread_pool_size() == 2
before = dict(os.environ)
try:
    vdjaqp.resolve_cpu_budget(1)
except RuntimeError as error:
    assert 'Polars' in str(error) and 'process' in str(error)
else:
    raise AssertionError('declared budget shrank while native pool stayed at two')
assert dict(os.environ) == before
assert vdjaqp.cpu_budget() == 2
assert vdjaqp.resolve_cpu_budget(2) == 2
'''
    result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, timeout=30)
    if result.returncode == 77:
        pytest.skip('needs two allocated cores to exercise late shrink')
    assert result.returncode == 0, result.stderr
