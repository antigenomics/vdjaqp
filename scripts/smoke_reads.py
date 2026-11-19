"""Run published InSilicoSeq once per tiny declared library geometry; validate gzip pairs."""
import argparse
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import sys
import time


def selected_modes(geometry='both'):
    """Resolve a declared library geometry without inspecting existing outputs."""
    if geometry not in ('both', 'amplicon', 'rnaseq'):
        raise ValueError('geometry must be both, amplicon or rnaseq')
    return ('amplicon', 'rnaseq') if geometry == 'both' else (geometry,)


def expected_pair_count(path):
    """InSilicoSeq's counts file specifies individual reads, not pairs."""
    total, seen = 0, set()
    for line in Path(path).read_text().splitlines():
        fields = line.split('\t')
        if len(fields) != 2 or not fields[0] or fields[0] in seen:
            raise ValueError('counts require unique template IDs and one integer read count')
        count = int(fields[1])
        if count < 0 or count % 2:
            raise ValueError('paired simulation requires nonnegative even individual-read counts')
        seen.add(fields[0])
        total += count // 2
    if total == 0:
        raise ValueError('at least one read pair is required')
    return total


def run(root, cpus, geometry='both', stage_timeout=300., *, pigz=None, seed=2468):
    if not math.isfinite(stage_timeout) or stage_timeout <= 0:
        raise ValueError('stage-timeout must be finite and positive')
    if type(seed) is not int or not 0 <= seed < 2**32:
        raise ValueError('read seed must be an integer in [0,2**32)')
    modes = selected_modes(geometry)
    root = Path(root)
    manifest = root/'read_run_manifest.json'
    if manifest.exists() or manifest.is_symlink():
        raise FileExistsError(f'read run manifest already exists: {manifest}')
    for scenario in modes:
        for suffix in ('R1', 'R2'):
            for extension in ('fastq', 'fastq.gz'):
                path = root/f'{scenario}_{suffix}.{extension}'
                if path.exists() or path.is_symlink():
                    raise FileExistsError(f'read output already exists: {path}')
    env=os.environ.copy()
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):env[name]='1'
    compression_identity = None
    if pigz is not None:
        if type(cpus) is not int or cpus < 1:
            raise ValueError('pigz requires a positive integer CPU budget')
        pigz = Path(pigz).resolve()
        if not pigz.is_file() or not os.access(pigz, os.X_OK):
            raise ValueError('--pigz must name an executable file')
        with pigz.open('rb') as stream:
            binary_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
        compression_identity = dict(binary=str(pigz), binary_sha256=binary_hash,
            version=subprocess.check_output([str(pigz), '-V'], stderr=subprocess.STDOUT,
                text=True, env=env, timeout=stage_timeout).strip())
    expected_pairs = expected_pair_count(root/'read_counts.tsv')
    records=[]
    binary=str(Path(sys.executable).parent/'iss')
    for scenario in modes:
        mode = 'amplicon' if scenario == 'amplicon' else 'metagenomics'
        prefix=root/scenario
        command=[binary,'generate','--genomes',str(root/f'{scenario}_templates.fasta'),'--readcount_file',str(root/'read_counts.tsv'),'--sequence_type',mode,'--model','MiSeq','--seed',str(seed),'--cpus',str(cpus),'--compress','--output',str(prefix)]
        if pigz is not None:command.remove('--compress')
        if scenario=='rnaseq':command+=['--fragment-length','400','--fragment-length-sd','40']
        start=time.perf_counter()
        with (root/f'{scenario}.log').open('w') as log:
            subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=stage_timeout)
        compression = None
        if pigz is not None:
            raw_paths = [root/f'{scenario}_{suffix}.fastq' for suffix in ('R1', 'R2')]
            raw_identity = {}
            for suffix, path in zip(('R1', 'R2'), raw_paths):
                if Path(str(path)+'.gz').exists() or Path(str(path)+'.gz').is_symlink():
                    raise FileExistsError(f'compressed read output already exists: {path}.gz')
                with path.open('rb') as stream:
                    digest = hashlib.file_digest(stream, 'sha256').hexdigest()
                raw_identity[suffix] = dict(bytes=path.stat().st_size, sha256=digest)
            # Keep both sources until both compressed mates pass decoded validation.
            compress_command = [str(pigz), '-k', '-n', '-p', str(cpus), '--', *map(str, raw_paths)]
            tick = time.perf_counter()
            with (root/f'{scenario}.log').open('a') as log:
                subprocess.run(compress_command, env=env, stdout=log, stderr=subprocess.STDOUT,
                               check=True, timeout=stage_timeout)
            compression = dict(**compression_identity, command=compress_command,
                               wall_seconds=time.perf_counter()-tick, logical_fastq=raw_identity)
        counts=[];logical_fastq={}
        for suffix in ('R1','R2'):
            path=root/f'{scenario}_{suffix}.fastq.gz'
            decoded_hash=hashlib.sha256();decoded_bytes=lines=0
            with gzip.open(path,'rb') as handle:
                for line in handle:
                    lines+=1
                    decoded_hash.update(line);decoded_bytes+=len(line)
            if compression is not None and dict(bytes=decoded_bytes,sha256=decoded_hash.hexdigest()) != raw_identity[suffix]:
                raise ValueError(f'compressed FASTQ differs from plaintext: {suffix}')
            assert lines%4==0
            counts.append(lines//4)
            logical_fastq[suffix]=dict(bytes=decoded_bytes,sha256=decoded_hash.hexdigest(),pairs=lines//4)
        assert counts==[expected_pairs,expected_pairs],counts
        record=dict(scenario=scenario,command=command,pairs=counts[0],wall_seconds=time.perf_counter()-start,peak_child_rss_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*(1 if sys.platform=='darwin' else 1024),model='MiSeq',cpus=cpus,seed=seed,logical_fastq=logical_fastq)
        if compression is not None:
            for path in raw_paths:path.unlink()
            record['compression']=dict(**compression,decoded_identity_verified=True,
                                       consumed_plaintext_removed=[p.name for p in raw_paths])
        records.append(record)
    (root/'read_run_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
    return records


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--cpus',type=int,required=True)
    p.add_argument('--geometry', choices=('both','amplicon','rnaseq'), default='both')
    p.add_argument('--stage-timeout',type=float,default=300.)
    p.add_argument('--read-seed',type=int,default=2468)
    p.add_argument('--pigz',type=Path,help='explicit native parallel gzip executable; default: ISS compression')
    a=p.parse_args()
    if not math.isfinite(a.stage_timeout) or a.stage_timeout<=0:p.error('stage-timeout must be finite and positive')
    print(json.dumps(run(a.root,a.cpus,a.geometry,a.stage_timeout,pigz=a.pigz,seed=a.read_seed),indent=2))
