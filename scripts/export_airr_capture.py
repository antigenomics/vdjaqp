"""Fresh ten-context agreeing-residue capture from supplied experimental inputs."""
import argparse,gc,hashlib,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from vdjaqp import resolve_cpu_budget


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--aggregation',choices=['alignment_union','exclusive_spectrum'],default='alignment_union')
    p.add_argument('--counts',type=Path);p.add_argument('--count-unit',choices=['reads','UMIs']);p.add_argument('--count-provenance')
    p.add_argument('--cpus',type=int,default=16);a=p.parse_args();resolve_cpu_budget(a.cpus)
    if a.output.exists():p.error('fresh output required')
    if sum(bool(x) for x in (a.counts,a.count_unit,a.count_provenance)) not in (0,3):p.error('counts need unit and provenance together')
    import polars as pl
    from vdjaqp.pipeline import Evidence
    from vdjaqp.data import load_targets,load_peptides,ms_stem_of
    from vdjaqp.airr_capture import capture_tables
    def sha(f):
        with f.open('rb') as s:return hashlib.file_digest(s,'sha256').hexdigest()
    refs=load_targets(a.input_root);single=refs.filter(pl.col('has_single_run')=='yes')
    if refs.height!=12 or single.height!=8:raise ValueError('expected twelve references/eight individual runs')
    specs=[('gel','dataset_mixa_gel',None,None),('solution','dataset_mixa',None,None)]+[(f'individual_{i+1:02d}','dataset_single_mabs',r['ms_stem'],r['clone_id']) for i,r in enumerate(single.iter_rows(named=True))]
    repo=Path(__file__).resolve().parents[1];source={str(f.relative_to(repo)):sha(f) for f in [*sorted((repo/'vdjaqp').glob('*.py')),Path(__file__).resolve()]}
    inputs={str(f.relative_to(a.input_root)):sha(f) for f in a.input_root.rglob('*') if f.is_file() and (f.suffix=='.parquet' or f.name=='target_clonotypes.csv')}
    count=pl.read_csv(a.counts,separator='\t') if a.counts else None
    a.output.mkdir(parents=True);report=dict(status='running',source_sha256=source,raw_sha256=inputs,aggregation=a.aggregation,counts_sha256=sha(a.counts) if a.counts else None,count_unit=a.count_unit,count_provenance=a.count_provenance,reference_counts='unmeasured; never zero imputation',cases=[])
    try:
        for case,dataset,stem,target in specs:
            observed=None
            if target is not None:
                observed=load_peptides(dataset,a.input_root);observed=observed.with_columns(ms_stem_of(observed)).filter(pl.col('ms_stem')==stem).drop('ms_stem')
            ev=Evidence.build(dataset,root=a.input_root,peptides=observed,exclude=() if target is None else {target});scores=ev.score(aggregation=a.aggregation)
            tables=capture_tables(ev,scores,set(refs['clone_id']),counts=count,count_unit=a.count_unit,count_provenance=a.count_provenance)
            dest=a.output/case;dest.mkdir()
            for name,frame in tables.items():frame.write_parquet(dest/(name+'.parquet'),compression='zstd')
            report['cases'].append(dict(case=case,records=ev.clonotypes.height,peptides=ev.peptides.height,alignments=ev.alignments.height))
            del tables,scores,ev,observed;gc.collect()
        if any(sha(repo/n)!=h for n,h in source.items()) or any(sha(a.input_root/n)!=h for n,h in inputs.items()):raise ValueError('source/input drift during export')
        report['status']='complete'
    except BaseException as e:report.update(status='failed',error=repr(e));raise
    finally:
        report['artifact_sha256']={str(f.relative_to(a.output)):sha(f) for f in a.output.rglob('*.parquet')};(a.output/'provenance.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
