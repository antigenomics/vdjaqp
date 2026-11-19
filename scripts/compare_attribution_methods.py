"""Fresh ten-context attribution comparison against both preserved final baselines.

Private local outputs only. No resume, computed feature reuse, threshold scan or
label-dependent scoring. Saved baseline scores are accepted requested deliverables.
"""
import argparse,gc,hashlib,json,subprocess,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from vdjaqp import resolve_cpu_budget


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root',type=Path,required=True)
    p.add_argument('--accepted-baselines',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cpus',type=int,default=16)
    a=p.parse_args();budget=resolve_cpu_budget(a.cpus)
    if a.output.exists():p.error('fresh output required')
    import numpy as np
    import polars as pl
    from vdjaqp.pipeline import Evidence
    from vdjaqp.bench import classification_metrics,rank_targets
    from vdjaqp.data import load_targets,load_peptides,load_clonotypes,ms_stem_of
    targets=load_targets(a.input_root);single=targets.filter(pl.col('has_single_run')=='yes')
    if targets.height!=12 or single.height!=8:p.error('expected twelve references/eight individual runs')
    specs=[('gel','dataset_mixa_gel',None,None),('solution','dataset_mixa',None,None)]
    specs += [(f'individual_{i+1:02d}','dataset_single_mabs',r['ms_stem'],r['clone_id']) for i,r in enumerate(single.iter_rows(named=True))]
    repo=Path(__file__).resolve().parents[1]
    sources={str(f.relative_to(repo)):sha(f) for f in [*sorted((repo/'vdjaqp').glob('*.py')),Path(__file__).resolve()]}
    inputs={str(f.relative_to(a.input_root)):sha(f) for f in sorted(a.input_root.rglob('*')) if f.is_file() and (f.suffix=='.parquet' or f.name=='target_clonotypes.csv')}
    baseline_files={f'{v}/{case}':a.accepted_baselines/v/case/'scores.parquet' for v in ['old','new'] for case,*_ in specs}
    baseline_hashes={k:sha(f) for k,f in baseline_files.items()}
    a.output.mkdir(parents=True)
    provenance=dict(status='running',cpu_budget=budget,source_sha256=sources,raw_sha256=inputs,
        accepted_final_score_sha256=baseline_hashes,baseline_labels={'old':'original arithmetic alignment_union','new':'corrected arithmetic alignment_union'},
        candidate='corrected arithmetic exclusive_spectrum',primary_metric='tie-group average_precision; not trapezoidal PR-AUC',
        threshold='confirmed_nats >= 20',positive_definition='12 reference records per mixture; one held-out target record per individual run',
        negative_definition='all other candidate records are benchmark negatives, not proven biological nonbinders',
        revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip())
    def save(): (a.output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    save();metrics=[];aliases=[]
    try:
        for case,dataset,stem,target in specs:
            start=time.monotonic();positives=set(targets['clone_id']) if target is None else {target}
            clones=load_clonotypes(dataset,a.input_root)
            if clones.height!=33052 or clones['clone_id'].n_unique()!=33052:raise ValueError('complete candidate universe changed')
            # Exact supplied-variable equality diagnostic, independent of benchmark labels.
            variable={}
            for row in clones.select('clone_id','seq','region_ids').iter_rows(named=True):
                seq=row['seq'];regions=row['region_ids'];n=sum(x==8 for x in regions)
                if n and regions[-n:]!=[8]*n:raise ValueError('nonterminal constant interval')
                variable[row['clone_id']]=(seq[:-n] if n else seq).upper().replace('I','L')
            target_classes={variable[k] for k in positives}
            for k,v in variable.items():
                if v in target_classes:aliases.append(dict(case=case,clone_id=k,record_positive=k in positives,full_supplied_variable_IL_alias=True))
            for variant in ['old','new']:
                scores=pl.read_parquet(baseline_files[f'{variant}/{case}'])
                if set(scores['clone_id'])!=set(clones['clone_id']):raise ValueError('baseline candidate membership differs')
                metrics.append(dict(case=case,method=provenance['baseline_labels'][variant],**classification_metrics(scores,positives)))
            observed=None
            if target is not None:
                observed=load_peptides(dataset,a.input_root)
                observed=observed.with_columns(ms_stem_of(observed)).filter(pl.col('ms_stem')==stem).drop('ms_stem')
            ev=Evidence.build(dataset,root=a.input_root,peptides=observed,exclude=() if target is None else {target})
            scores=ev.score(aggregation='exclusive_spectrum')
            if set(scores['clone_id'])!=set(clones['clone_id']):raise ValueError('variant candidate membership differs')
            dest=a.output/case;dest.mkdir();scores.write_parquet(dest/'scores.parquet',compression='zstd')
            rank_targets(scores,'score',positives).write_parquet(dest/'target_ranks.parquet',compression='zstd')
            metrics.append(dict(case=case,method=provenance['candidate'],**classification_metrics(scores,positives)))
            pl.DataFrame(metrics).write_csv(a.output/'metrics.tsv',separator='\t')
            print(json.dumps(dict(case=case,seconds=time.monotonic()-start,peptides=ev.peptides.height,alignments=ev.alignments.height)),flush=True)
            del ev,scores,observed,clones;gc.collect()
        pl.DataFrame(aliases).write_csv(a.output/'private_target_sequence_aliases.tsv',separator='\t')
        if any(sha(repo/n)!=h for n,h in sources.items()):raise ValueError('source changed during comparison')
        if any(sha(a.input_root/n)!=h for n,h in inputs.items()):raise ValueError('raw inputs changed during comparison')
        if any(sha(baseline_files[n])!=h for n,h in baseline_hashes.items()):raise ValueError('accepted baseline changed')
        provenance['status']='complete'
    except BaseException as error:
        provenance.update(status='failed',error=repr(error));raise
    finally:
        provenance['artifact_sha256']={str(f.relative_to(a.output)):sha(f) for f in a.output.rglob('*') if f.is_file() and f.name!='provenance.json'};save()


if __name__=='__main__':main()
