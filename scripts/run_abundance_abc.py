"""Run pyABC only after one qualified fresh joint draw; raw inputs stay local.

The scientific forward worker always uses the baseline .venv, not inference-env
versions of NumPy/ARDA/scoring dependencies. No previous summary is a model cache.
"""
import argparse,csv,hashlib,json,os,subprocess,sys,time
from pathlib import Path
repo=Path(__file__).resolve().parents[1];sys.path.insert(0,str(repo))
from vdjaqp import resolve_cpu_budget


def verify_sources(qualified,source_root,reconciliation_path,qualification_report):
    """Permit only explicit reviewed source pairs, each pinned to an approved snapshot."""
    def sha(path):
        with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
    reconciliation_path=Path(reconciliation_path)
    record=json.loads(reconciliation_path.read_text())
    allowed={'vdjaqp/abundance_inference.py','vdjaqp/transcript_boundary.py','vdjaqp/transcript_context.py'}
    exceptions=record['source_exceptions']
    if set(exceptions)!=allowed or record['qualification_report_sha256']!=sha(qualification_report):
        raise ValueError('unrecognised source exception or qualification report')
    receipt=dict(reconciliation_sha256=sha(reconciliation_path),qualification_report_sha256=sha(qualification_report),
                 source_pairs={},approved_snapshot_sha256={},verified_current_sources={})
    for name,expected in qualified['source_sha256'].items():
        actual=sha(source_root/name)
        if name in exceptions:
            item=exceptions[name]
            if item['qualified_sha256']!=expected or item['approved_current_sha256']!=actual:
                raise ValueError('reviewed source hash pair no longer matches: '+name)
            snapshot=Path(item['approved_snapshot'])
            if snapshot.is_absolute() or '..' in snapshot.parts:
                raise ValueError('approved snapshot must be local to reconciliation')
            if sha(reconciliation_path.parent/snapshot)!=actual:
                raise ValueError('approved snapshot hash mismatch: '+name)
            receipt['source_pairs'][name]=dict(qualified_sha256=expected,approved_current_sha256=actual)
            receipt['approved_snapshot_sha256'][str(snapshot)]=actual
            if item.get('qualified_snapshot'):
                old=Path(item['qualified_snapshot'])
                if old.is_absolute() or '..' in old.parts or sha(reconciliation_path.parent/old)!=expected:
                    raise ValueError('qualified source snapshot hash mismatch')
                receipt['approved_snapshot_sha256'][str(old)]=expected
        elif actual!=expected:
            raise ValueError('qualified scientific dependency changed: '+name)
        receipt['verified_current_sources'][name]=actual
    if not allowed<=receipt['source_pairs'].keys():raise ValueError('reviewed source absent from qualification')
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repertoire',type=Path,required=True);p.add_argument('--target-ids',type=Path,required=True)
    p.add_argument('--targets',type=Path,required=True);p.add_argument('--source-format',required=True)
    p.add_argument('--qualified-draw',type=Path,required=True);p.add_argument('--qualification-resource',type=Path,required=True)
    p.add_argument('--source-reconciliation',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--mssimulator',type=Path,required=True)
    p.add_argument('--cpus',type=int,default=16)
    a=p.parse_args();budget=resolve_cpu_budget(a.cpus)
    if a.output.exists():p.error('fresh output required; no inference resume')
    os.environ['MPLCONFIGDIR']=str(repo/'.tool-envs/inference/mpl')
    os.environ['XDG_CACHE_HOME']=str(repo/'.tool-envs/inference/cache')
    from vdjaqp.abundance_inference import target_vector,run_abc,distance
    rows=list(csv.DictReader(a.targets.open(),delimiter='\t'))
    for row in rows:
        for field in ('value','numerator','denominator'):
            row[field]=None if row[field] in ('','null') else float(row[field])
    observed=target_vector(rows)
    resource=json.loads(a.qualification_resource.read_text())
    qualified=json.loads((a.qualified_draw/'joint_draw.json').read_text())
    if resource['failure'] is not None or resource['exit_code']!=0 or qualified['status']!='complete':
        p.error('first joint draw must pass resource/producer qualification')
    def sha(path):
        with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
    for key,path in [('repertoire',a.repertoire),('target_ids',a.target_ids)]:
        if qualified['raw_input_sha256'][key]!=sha(path):p.error('qualified raw inputs differ')
    if len(qualified['relative_amounts'])!=12 or any(abs(v-1)>1e-12 for v in qualified['relative_amounts']):
        p.error('initial qualification must be the declared equal-input engineering draw')
    if any(r['selected_mssimulator_sha256']!=sha(a.mssimulator) for r in qualified['native_reports'].values()):
        p.error('native executable differs from qualified draw')
    try:
        verified_sources=verify_sources(qualified,repo,a.source_reconciliation,a.qualified_draw/'joint_draw.json')
    except ValueError as error:p.error(str(error))
    if distance(json.loads((a.qualified_draw/'summary.json').read_text()),observed)<0:p.error('invalid qualified distance')
    remaining=21600.-resource['wall_seconds']
    if remaining<=0:p.error('qualification consumed aggregate six-hour starting bound')
    draws=a.output/'draws'
    def forward(amounts,index,seconds):
        draws.mkdir(parents=True,exist_ok=True)
        amount_file=draws/f'amounts-{index:04d}.json';amount_file.write_text(json.dumps(amounts.tolist())+'\n')
        case=draws/f'draw-{index:04d}'
        command=[str(repo/'.venv/bin/python'),str(repo/'scripts/profile_build.py'),
            '--output',str(draws/f'resource-{index:04d}.json'),'--project-root',str(repo),'--workload-root',str(a.output.resolve()),
            '--timeout',str(min(1800.,seconds)),'--rss-limit','24000000000','--disk-limit','10000000000',
            '--cpus',str(budget),'--',str(repo/'.venv/bin/python'),str(repo/'scripts/run_abundance_forward.py'),
            '--repertoire',str(a.repertoire.resolve()),'--target-ids',str(a.target_ids.resolve()),'--targets',str(a.targets.resolve()),
            '--source-format',a.source_format,'--output',str(case.resolve()),'--amounts',str(amount_file.resolve()),
            '--cpus',str(budget),'--mssimulator',str(a.mssimulator.resolve())]
        subprocess.run(command,cwd=repo,check=True,timeout=min(1800.,seconds)+60)
        return json.loads((case/'summary.json').read_text())
    try:
        report=run_abc(forward,observed,a.output,max_draws=255,wall_seconds=remaining)
    finally:
        if a.output.is_dir():
            (a.output/'run_provenance.json').write_text(json.dumps(dict(
                raw_input_sha256={key:sha(path) for key,path in [('repertoire',a.repertoire),('target_ids',a.target_ids),('targets',a.targets)]},
                qualification_report_sha256=sha(a.qualified_draw/'joint_draw.json'),
                verified_source_reconciliation=verified_sources,
                qualification_resource_sha256=sha(a.qualification_resource),
                code_sha256={name:sha(repo/name) for name in ['scripts/run_abundance_abc.py','scripts/run_abundance_forward.py','vdjaqp/abundance_inference.py','vdjaqp/abundance_forward.py']},
                inference_lock_sha256=sha(repo/'envs/inference/requirements.lock'),
                first_qualification_count=1,total_attempt_bound=256,fresh_attempt_bound=255,
                qualification_wall_seconds=resource['wall_seconds'],fresh_wall_bound_seconds=remaining,
                input_intensity=10000,ms2_per_bin=5,source_format=a.source_format,
                sigma='ln(10)',coordinates=255,recovery_used=False),indent=2)+'\n')
    if report['status']!='complete':raise SystemExit(1)

if __name__=='__main__':main()
