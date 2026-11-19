"""Published fastp overlap merging only; preserve unmerged input mates."""
import gzip
import json
from pathlib import Path
import subprocess
import time


def preprocess_pairs(r1, r2, output, *, threads, fastp, stage_timeout=None):
    output=Path(output)
    if output.exists():raise FileExistsError(f'preprocessing output exists: {output}')
    report=output.with_suffix('.json');html=output.with_suffix('.html')
    command=[str(fastp),'--in1',str(r1),'--in2',str(r2),'--merge','--include_unmerged',
        '--merged_out',str(output),'--json',str(report),'--html',str(html),
        '--thread',str(threads),'--disable_adapter_trimming','--disable_trim_poly_g',
        '--disable_quality_filtering','--disable_length_filtering']
    start=time.perf_counter();subprocess.run(command,check=True,timeout=stage_timeout)
    stats=json.loads(report.read_text())
    def count(path):
        with gzip.open(path,'rt') as f:
            lines=sum(1 for _ in f)
        if lines%4:raise ValueError('invalid FASTQ record count')
        return lines//4
    first,second=count(r1),count(r2)
    if first!=second:raise ValueError('input paired read counts differ')
    records=count(output)
    with gzip.open(output,'rt') as f:
        merged=sum(1 for i,line in enumerate(f) if i%4==0 and any(x.startswith('merged_') for x in line.split()[1:]))
    if stats['filtering_result']['passed_filter_reads']!=2*first:
        raise ValueError('native preprocessing unexpectedly filtered input reads')
    if stats['merged_and_filtered']['total_reads']!=records:
        raise ValueError('native output count disagrees with FASTQ')
    if records!=merged+2*(first-merged):raise ValueError('merged/unmerged retention accounting differs')
    return {'command':command,'wall_seconds':time.perf_counter()-start,
        'input_pairs':first,'merged_pairs':merged,'unmerged_pairs':first-merged,
        'retained_output_records':records,'filters_disabled':True}
