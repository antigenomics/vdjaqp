"""One fixed CPU budget per run, resolved before any numerical library loads."""

from __future__ import annotations

import os
import sys
from numbers import Integral

__all__ = ["DEFAULT_CPU_BUDGET", "resolve_cpu_budget", "cpu_budget", "bounded_threads"]

#: Use the actual available allocation by default.
DEFAULT_CPU_BUDGET = "all"

#: Compute engines that own the allocated parallelism.
_THREAD_VARS = (
    "POLARS_MAX_THREADS", "RAYON_NUM_THREADS", "OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT", "NUMEXPR_NUM_THREADS",
)

_RESOLVED: int | None = None


def _available_cores() -> int:
    """Affinity/cgroup ceiling from vdjtools, further limited by SLURM."""
    from vdjtools.cores import available_cores

    available = available_cores()
    allocation = os.environ.get("SLURM_CPUS_PER_TASK")
    if allocation is not None:
        try:
            requested = int(allocation)
        except ValueError:
            raise ValueError("SLURM_CPUS_PER_TASK must be a positive integer") from None
        if requested < 1:
            raise ValueError("SLURM_CPUS_PER_TASK must be a positive integer")
        available = min(available, requested)
    return available


def resolve_cpu_budget(budget: int | str | None = None) -> int:
    """Fix the ceiling for native compute engines; keep auxiliary BLAS serial.

    Args:
        budget: ``None`` for :data:`DEFAULT_CPU_BUDGET`, ``"serial"`` or ``1``
            for a single compute thread, ``"all"`` for every allocated core, or
            a positive integer.

    Returns:
        The resolved budget ``B``.

    Raises:
        ValueError: On a non-positive or unparseable request. Positive requests
            above the available allocation are clamped to that allocation.
        RuntimeError: If an already initialized Polars pool exceeds the new
            ceiling. No environment or declared budget is changed on failure.
    """
    global _RESOLVED
    available = _available_cores()
    if budget is None:
        resolved = available
    elif budget == "serial":
        resolved = 1
    elif budget == "all":
        resolved = available
    else:
        if isinstance(budget, bool) or not isinstance(budget, (Integral, str)):
            raise ValueError("cpu budget must be an integer, serial or all")
        try:
            resolved = int(budget)
        except (TypeError, ValueError):
            raise ValueError(
                f"cpu budget must be an int, 'serial' or 'all', got {budget!r}"
            ) from None
        if resolved < 1:
            raise ValueError(f"cpu budget must be >= 1, got {resolved}")
        resolved = min(resolved, available)
    previous = {name: os.environ.get(name) for name in _THREAD_VARS}
    for name in _THREAD_VARS:
        os.environ[name] = str(resolved)
    try:
        polars = sys.modules.get("polars")
        if polars is not None and polars.thread_pool_size() > resolved:
            raise RuntimeError(
                "Polars already has more threads than the requested budget; start a fresh "
                "process and resolve its CPU budget before importing Polars"
            )
    except BaseException:
        for name, value in previous.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value
        raise
    # Keep the small 20-letter background mixing product serial; search/table
    # engines retain the full budget. Avoid idle BLAS teams during native stages.
    for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[name] = "1"
    # Polars' separate Tokio executor otherwise adds up to 32 idle workers.
    # Keep its I/O pool bounded while retaining the full compute allocation.
    os.environ["POLARS_ASYNC_THREAD_COUNT"] = "1"
    _RESOLVED = resolved
    return resolved


def cpu_budget() -> int:
    """The budget resolved for this process."""
    if _RESOLVED is None:
        return resolve_cpu_budget()
    return _RESOLVED


def bounded_threads(threads: int | None = None) -> int:
    """Per-call native threads, capped by the process budget without changing it.

    ``None`` and zero use the selected process budget. Positive integers may
    request less; oversized requests are clamped. Invalid types/negatives fail.
    """
    if threads is None:
        return cpu_budget()
    if isinstance(threads, bool) or not isinstance(threads, Integral) or threads < 0:
        raise ValueError("threads must be a non-negative integer or None")
    return cpu_budget() if threads == 0 else min(int(threads), cpu_budget())
