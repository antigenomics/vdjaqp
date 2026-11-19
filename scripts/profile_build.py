"""Measure and bound a complete local synthetic build's process tree.

Reports sampled sums of RSS (not PSS), live threads, disk footprint and wall time.
Only the newly launched process group can be terminated by this wrapper.
"""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import traceback


def process_scope(rows, root_pid):
    """Select the owned process group and descendants, each PID once."""
    selected = {root_pid} | {r[0] for r in rows if r[2] == root_pid}
    while True:
        expanded = selected | {r[0] for r in rows if r[1] in selected}
        if expanded == selected:
            return [r for r in rows if r[0] in selected]
        selected = expanded


def snapshot(root_pid):
    result = subprocess.run(['ps','-axo','pid=,ppid=,pgid=,rss=,pcpu='],
                            capture_output=True,text=True,check=True,timeout=10)
    rows=[]
    for line in result.stdout.splitlines():
        pid,parent,group,rss,cpu=line.split()
        rows.append((int(pid),int(parent),int(group),int(rss),float(cpu)))
    owned=process_scope(rows,root_pid)
    pids=','.join(str(r[0]) for r in owned if r[3]>0)
    threads=0
    if pids:
        command=['ps','-M','-p',pids] if sys.platform=='darwin' else ['ps','-L','-p',pids,'-o','pid=']
        listed=subprocess.run(command,capture_output=True,text=True,timeout=10)
        if listed.returncode not in (0,1):
            raise RuntimeError('native thread enumeration failed')
        threads=max(0,len(listed.stdout.splitlines())-(1 if sys.platform=='darwin' else 0))
    return owned,threads


def terminate_owned(proc):
    """Stop this launch's process group, never the caller's group."""
    try:os.killpg(proc.pid,signal.SIGTERM)
    except ProcessLookupError:pass
    try:proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:os.killpg(proc.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        proc.wait(timeout=10)


def du_allocated(path):
    """Allocated bytes and whether vanished descendants made this scan partial."""
    path = Path(path).resolve()
    if not path.exists():
        return 0, False
    result = subprocess.run(['du', '-sk', str(path)], capture_output=True, text=True,
                            timeout=30, env=dict(os.environ, LC_ALL='C'))
    lines = result.stdout.splitlines()
    fields = lines[0].split(maxsplit=1) if len(lines) == 1 else []
    valid = len(fields) == 2 and fields[0].isascii() and fields[0].isdigit() and fields[1] == str(path)
    partial = result.returncode == 1 and bool(result.stderr.strip()) and path.is_dir()
    for line in result.stderr.splitlines():
        prefix, suffix = 'du: ', ': No such file or directory'
        if not line.startswith(prefix) or not line.endswith(suffix):
            partial = False
            break
        child = line[len(prefix):-len(suffix)]
        for operation in ('cannot access ', 'cannot read directory '):
            if child.startswith(operation):
                child = child[len(operation):]
                break
        child = Path(child.strip("'"))
        if not child.is_absolute() or child == path or not child.is_relative_to(path):
            partial = False
            break
    if not valid or (result.returncode != 0 and not partial) or (result.returncode == 0 and result.stderr):
        raise RuntimeError(f'du measurement failed: returncode={result.returncode}; '
                           f'stdout={result.stdout!r}; stderr={result.stderr!r}')
    return int(fields[0]) * 1024, partial


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--project-root',type=Path,default=Path.cwd())
    p.add_argument('--workload-root',type=Path,required=True,help='sole directory receiving generated job outputs')
    p.add_argument('--timeout',type=float,default=1800)
    p.add_argument('--rss-limit',type=int,default=24_000_000_000)
    p.add_argument('--disk-limit',type=int,default=10_000_000_000)
    p.add_argument('--cpus',type=int)
    p.add_argument('command',nargs=argparse.REMAINDER)
    a=p.parse_args()
    command=a.command[1:] if a.command[:1]==['--'] else a.command
    logpath=a.output.with_suffix('.log')
    if not command or a.timeout<=0 or a.rss_limit<=0 or a.disk_limit<=0:
        p.error('a command and positive resource bounds are required')
    if a.output.exists() or logpath.exists():
        p.error('report/log already exists; select a new output')
    if not a.project_root.is_dir():
        p.error('project root must exist')
    if not a.workload_root.resolve().is_relative_to(a.project_root.resolve()):
        p.error('workload root must be within the measured project')
    repo=Path(__file__).resolve().parents[1]
    sys.path.insert(0,str(repo))
    from vdjaqp import resolve_cpu_budget
    budget=resolve_cpu_budget(a.cpus)
    # The parent retains idle Polars and native Rayon pools during sequential
    # external stages. Allow those plus one active child pool and helper threads.
    thread_limit=3*budget+16
    # Verify monitoring access before starting any child work.
    a.output.parent.mkdir(parents=True,exist_ok=True)
    partial_scans = 0
    def allocated(path):
        nonlocal partial_scans
        size, partial = du_allocated(path)
        partial_scans += int(partial)
        return size
    try:
        snapshot(os.getpid())
        baseline_disk=allocated(a.project_root)
        initial_workload=allocated(a.workload_root)
    except BaseException as error:
        report = dict(command=command, cpu_budget=budget, exit_code=None,
                      failure='preflight_monitor_error', samples=0,
                      peak_sampled_summed_rss_bytes=0, peak_sampled_project_allocated_bytes=None,
                      transient_partial_disk_scans=partial_scans,
                      monitor_error=dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc()))
        a.output.write_text(json.dumps(report, indent=2)+'\n')
        raise SystemExit(1) from error
    if baseline_disk>a.disk_limit:p.error('project already exceeds disk limit')
    peak_rss=peak_threads=peak_processes=peak_disk=0
    peak_cpu=0.
    samples=0
    failure=None
    monitor_error=None
    start=time.monotonic()
    previous_sample=start
    max_gap=0.
    with logpath.open('w') as log:
        proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,
                              start_new_session=True,cwd=a.project_root)
        try:
            while True:
                rows,threads=snapshot(proc.pid)
                now=time.monotonic()
                max_gap=max(max_gap,now-previous_sample)
                previous_sample=now
                samples+=1
                rss=sum(r[3] for r in rows)*1024
                peak_rss=max(peak_rss,rss)
                peak_threads=max(peak_threads,threads)
                peak_processes=max(peak_processes,len(rows))
                peak_cpu=max(peak_cpu,sum(r[4] for r in rows))
                elapsed=time.monotonic()-start
                # Immutable runtimes dominate project files; sample only this
                # job's output growth, then recheck the whole project at exit.
                log_bytes=logpath.stat().st_blocks*512
                peak_disk=max(peak_disk,baseline_disk-initial_workload+allocated(a.workload_root)+log_bytes)
                if rss>a.rss_limit:failure='summed_rss_limit'
                elif peak_disk>a.disk_limit:failure='project_disk_limit'
                elif threads>thread_limit:failure='live_thread_limit'
                elif elapsed>a.timeout:failure='whole_build_timeout'
                if failure:
                    log.flush()
                    with logpath.open('a') as diagnostic:
                        diagnostic.write('\nResource-stop process rows (pid,parent,group,RSS KiB,CPU %):\n'+repr(rows)+'\n')
                        if sys.platform=='darwin' and rows:
                            subprocess.run(['ps','-M','-p',','.join(str(r[0]) for r in rows)],
                                           stdout=diagnostic,stderr=diagnostic,timeout=10)
                    break
                if proc.poll() is not None:
                    remaining,_=snapshot(proc.pid)
                    if any(r[3]>0 for r in remaining):failure='live_children_after_command_exit'
                    break
                try:proc.wait(timeout=.25)
                except subprocess.TimeoutExpired:pass
        except BaseException as error:
            failure = 'monitor_error'
            monitor_error = dict(type=type(error).__name__, message=str(error),
                                 traceback=traceback.format_exc())
        if failure:
            try:
                terminate_owned(proc)
            except BaseException as error:
                termination_error = dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc())
                if monitor_error is None:
                    failure = 'monitor_error'
                    monitor_error = termination_error
                else:
                    monitor_error['termination_error'] = termination_error
        status = proc.poll()
        if status is None and not failure:
            status = proc.wait()
    final_disk = None
    try:
        final_disk = allocated(a.project_root)
    except BaseException as error:
        if monitor_error is None:
            failure = 'monitor_error'
            monitor_error = dict(type=type(error).__name__, message=str(error),
                                 traceback=traceback.format_exc())
        else:
            monitor_error['final_disk_error'] = dict(type=type(error).__name__, message=str(error))
    if final_disk is not None and final_disk>a.disk_limit and failure is None:failure='project_disk_limit'
    report=dict(command=command,cpu_budget=budget,exit_code=status,failure=failure,
                owned_process_group=proc.pid,
                wall_seconds=time.monotonic()-start,samples=samples,sampling_pause_seconds=.25,
                maximum_observed_sample_gap_seconds=max_gap,
                peak_sampled_summed_rss_bytes=peak_rss,peak_live_threads=peak_threads,
                peak_process_count=peak_processes,peak_reported_cpu_percent=peak_cpu,
                peak_sampled_project_allocated_bytes=max(peak_disk,final_disk or 0),
                final_project_allocated_bytes=final_disk,
                monitor_error=monitor_error,transient_partial_disk_scans=partial_scans,
                partial_disk_scan_meaning='du ENOENT for vanished descendants only; root remained present and numeric summary retained; partial samples can underestimate concurrent disk use',
                disk_metric='preflight whole project plus declared workload/log growth; whole project rechecked at exit',
                memory_metric='sampled sum of process RSS; shared pages counted per process; not PSS',
                cpu_metric='sum of ps lifetime CPU percentages; not instantaneous utilization',
                live_thread_limit=thread_limit)
    a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if status or failure:raise SystemExit(status or 1)


if __name__=='__main__':main()
