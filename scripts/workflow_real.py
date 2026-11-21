"""Fresh full-candidate ten-context evaluation using existing Evidence APIs.

Truth is used only for declared calibration/exclusion and post-score evaluation.
Original producer records remain unchanged; no native ARDA result is claimed.
"""
import argparse
import gc
import hashlib
import json
import sys
from pathlib import Path
import time

from vdjaqp.data_hub import sha256 as sha


def export_mappings(evidence, output, query_catalog, *, query_role='identity_eligible'):
    """Retain seqtree placements independently of evidence/aggregation outputs."""
    import importlib.metadata
    output=Path(output)
    output.mkdir(parents=True,exist_ok=False)
    native=evidence.matches.with_row_index('placement_id')
    native.select('placement_id','clone_id','peptide_id','start','end','n_subs').write_parquet(output/'placements.parquet',compression='zstd')
    native.select('placement_id','mismatch_positions').explode('mismatch_positions',empty_as_null=True).drop_nulls('mismatch_positions').rename({'mismatch_positions':'peptide_offset'}).write_parquet(output/'substitutions.parquet',compression='zstd')
    evidence.clonotypes.select('clone_id','seq').write_parquet(output/'queries.parquet',compression='zstd')
    evidence.peptides.select('peptide_id','seq').write_parquet(output/'observations.parquet',compression='zstd')
    receipt=dict(status='complete',query_role=query_role,engine='seqtree.TextIndex via existing vdjaqp.match.PeptideMatcher',seqtree_version=importlib.metadata.version('seqtree'),query_catalog=str(query_catalog),query_catalog_sha256=sha(query_catalog),query_records=evidence.clonotypes.height,observations=evidence.peptides.height,placements=evidence.matches.height,contains_evidence_or_scores=False,coordinates='zero-based half-open in original query protein; substitution offsets relative to peptide; I/L mass equivalence preserved',artifact_sha256={path.name:sha(path) for path in sorted(output.glob('*.parquet'))})
    (output/'mapping_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def export_partial_mappings(query_catalog, peptides, output, cpus):
    """Map observed partial intervals without constructing evidence or scores."""
    from types import SimpleNamespace
    import polars as pl
    from vdjaqp.match import PeptideMatcher, shipped_budget
    queries=pl.read_parquet(query_catalog)
    eligible=queries['identity_eligible']
    if eligible.dtype != pl.Boolean or eligible.null_count() or eligible.any():
        raise ValueError('compatibility-only queries require nonmissing false identity_eligible')
    matcher=PeptideMatcher.from_clonotypes(queries)
    matches=matcher.search(peptides['peptide_id'],peptides['seq'],budget=shipped_budget,threads=cpus)
    return export_mappings(SimpleNamespace(clonotypes=queries,peptides=peptides,matches=matches),
                           output,query_catalog,query_role='partial_compatibility_only')


def build_evidence(dataset, input_root, calibration_root, peptides, target):
    """Keep query availability separate from an explicitly shared calibration."""
    from vdjaqp.pipeline import Evidence
    excluded=() if target is None else {target}
    if calibration_root is None:
        return Evidence.build(dataset,root=input_root,peptides=peptides,exclude=excluded)
    from vdjaqp.data import load_clonotypes,load_peptides
    from vdjaqp.errormodel import ErrorModel
    model=ErrorModel.calibrate(root=calibration_root,exclude=set(excluded))
    if peptides is None:peptides=load_peptides(dataset,input_root)
    return Evidence.from_frames(load_clonotypes(dataset,input_root),peptides,
                                error_model=model,dataset=dataset)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--mapping-output',type=Path,help='Separate fresh directory for native peptide placements, before score aggregation')
    p.add_argument('--calibration-root',type=Path,help='Explicit common calibration population; default uses each input root as before')
    p.add_argument('--compatibility-queries',type=Path,help='Observed partial query catalog: separate placements only, never scores')
    p.add_argument('--cpus',type=int)
    p.add_argument('--placeholder-count-one',action='store_true')
    p.add_argument('--expected-candidates',type=int,default=33052,help='Exact declared candidate population; default preserves the accepted cohort gate')
    p.add_argument('--case',choices=['gel','solution']+[f'individual_{i:02d}' for i in range(1,9)],help='One scheduled context; aggregate all ten before cohort acceptance')
    p.add_argument('--publication',action='store_true',help='also export the five existing method layers and complete peptide support')
    a=p.parse_args()
    if a.expected_candidates<=0:p.error('expected candidate count must be positive')
    import vdjaqp
    budget=vdjaqp.resolve_cpu_budget(a.cpus)
    import polars as pl
    from vdjaqp.data import load_targets,load_peptides,ms_stem_of
    from vdjaqp.airr_capture import capture_tables
    from vdjaqp.report import peptide_support_table
    from vdjaqp.data import isotype_of
    from vdjaqp.aggregate import isotype_log_prior
    from vdjaqp.aggregate import ClonotypeScorer
    from vdjaqp.bench import classification_metrics,rank_targets,summarise
    if a.output.exists():p.error('fresh output required')
    if a.mapping_output is not None and a.mapping_output.exists():p.error('fresh mapping output required')
    if a.compatibility_queries is not None and a.mapping_output is None:p.error('partial mappings require a separate mapping output')
    a.output.mkdir(parents=True)
    targets=load_targets(a.input_root)
    single=targets.filter(pl.col('has_single_run')=='yes')
    if targets.height!=12 or single.height!=8:raise ValueError('expected 12 references/8 individual runs')
    cases=[('gel','dataset_mixa_gel',None,None),('solution','dataset_mixa',None,None)]
    cases += [(f'individual_{i+1:02d}','dataset_single_mabs',r['ms_stem'],r['clone_id']) for i,r in enumerate(single.iter_rows(named=True))]
    if a.case is not None:cases=[case for case in cases if case[0]==a.case]
    metrics=[];layers=[]
    repo=Path(__file__).resolve().parents[1]
    source={str(p.relative_to(repo)):sha(p) for p in [*sorted((repo/'vdjaqp').glob('*.py')),Path(__file__).resolve()]}
    raw={str(p.relative_to(a.input_root)):sha(p) for p in sorted(a.input_root.rglob('*')) if p.is_file() and (p.suffix=='.parquet' or p.name=='target_clonotypes.csv')}
    begun=time.perf_counter()
    calibration_raw={} if a.calibration_root is None else {str(p.relative_to(a.calibration_root)):sha(p) for p in sorted(a.calibration_root.rglob('*')) if p.is_file() and (p.suffix=='.parquet' or p.name=='target_clonotypes.csv')}
    partial_sha=None if a.compatibility_queries is None else sha(a.compatibility_queries)
    report=dict(compatibility_query_sha256=partial_sha,calibration_root=str(a.calibration_root) if a.calibration_root is not None else None,calibration_raw_sha256=calibration_raw,status='running',source_sha256=source,raw_sha256=raw,cpus=budget,command=sys.argv,placeholder_counts=a.placeholder_count_one,count_unit='placeholder' if a.placeholder_count_one else None,cases=[])
    (a.output/'provenance.json').write_text(json.dumps(report,indent=2)+'\n')
    for case,dataset,stem,target in cases:
        start=time.perf_counter()
        observed=None
        if target:
            observed=load_peptides(dataset,a.input_root)
            observed=observed.with_columns(ms_stem_of(observed)).filter(pl.col('ms_stem')==stem).drop('ms_stem','source_file')
        ev=build_evidence(dataset,a.input_root,a.calibration_root,observed,target)
        model_parameters=dict(conf_x=ev.error._conf_x.tolist(),conf_e=ev.error._conf_e.tolist(),substitution=ev.error.substitution.tolist(),n_substitutions=ev.error.n_substitutions,e_db=ev.error.e_db)
        model_sha256=hashlib.sha256(json.dumps(model_parameters,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        if ev.clonotypes.height!=a.expected_candidates:raise ValueError('full declared candidate population required')
        mapping=None
        if a.mapping_output is not None:
            mapping_path=a.mapping_output if a.case is not None else a.mapping_output/case
            mapping=export_mappings(ev,mapping_path,a.input_root/dataset/'clonotypes.parquet')
        positives=set(targets['clone_id']) if target is None else {target}
        dest=a.output/case;dest.mkdir()
        counts=ev.clonotypes.select('clone_id').with_columns(count=pl.lit(1,dtype=pl.UInt64)) if a.placeholder_count_one else None
        # One current-run scorer owns the common residue information. Method
        # layers change aggregation options, not the candidate/background model.
        scorer=ClonotypeScorer(ev.clonotypes,ev.background)
        full_scores=None
        if a.publication:
            ev.alignments.write_parquet(dest/'alignment_evidence.parquet',compression='zstd')
            ev.clonotypes.with_columns(count=pl.lit(1,dtype=pl.UInt64),count_is_placeholder=pl.lit(True)).write_parquet(dest/'candidates.parquet',compression='zstd') if a.placeholder_count_one else ev.clonotypes.write_parquet(dest/'candidates.parquet',compression='zstd')
            peptide_support_table(ev).write_parquet(dest/'peptide_support.parquet',compression='zstd')
            baseline=ev.match_counts()
            isotypes=ev.clonotypes.select(isotype_of(ev.clonotypes))['isotype']
            isotype=ev.clonotypes.select('clone_id').with_columns(pl.Series('isotype_baseline',isotype_log_prior(isotypes.to_numpy())*1e6 if 'unknown' in isotypes else (isotypes=='IgG').cast(pl.Float64)*1e6))
            full_scores=scorer.score(ev.alignments,ev.peptides)
            methods=[('match_count',baseline),('match_count_isotype',baseline.join(isotype,on='clone_id').with_columns(score=pl.col('score')+pl.col('isotype_baseline'))),('no_competition',scorer.score(ev.alignments,ev.peptides,use_prior=False,competition=False)),('competition',scorer.score(ev.alignments,ev.peptides,use_prior=False)),('full',full_scores)]
            for name,table in methods:
                ranks=rank_targets(table,'score',positives)
                ranks.write_csv(dest/(name+'_target_ranks.tsv'),separator='\t')
                qualified=table.filter(pl.col('confirmed_nats')>=20) if 'confirmed_nats' in table.columns else None
                layers.append(dict(case=case,method=name,n_targets=len(positives),**summarise(ranks['rank'].to_numpy(),table.height),shortlist_count=qualified.height if qualified is not None else None,shortlist_target_count=qualified.filter(pl.col('clone_id').is_in(positives)).height if qualified is not None else None,framework_only_count=int(table['framework_only'].sum()) if 'framework_only' in table.columns else None,framework_only_target_count=int(table.filter(pl.col('clone_id').is_in(positives))['framework_only'].sum()) if 'framework_only' in table.columns else None,n_peptides=ev.peptides.height,n_alignments=ev.matches.height))
            pl.DataFrame(layers,infer_schema_length=None).write_csv(a.output/'method_layers.tsv',separator='\t')
            del methods,baseline
        for method in ['alignment_union','exclusive_spectrum']:
            scores=full_scores if method=='alignment_union' and full_scores is not None else scorer.score(ev.alignments,ev.peptides,aggregation=method)
            scores.write_parquet(dest/(method+'_scores.parquet'),compression='zstd')
            rank_targets(scores,'score',positives).write_csv(dest/(method+'_ranks.tsv'),separator='\t')
            if method=='alignment_union':
                for name,frame in capture_tables(ev,scores,set(targets['clone_id']),counts=counts,count_unit='placeholder' if counts is not None else None,count_provenance='author-requested count=1; not measured reads' if counts is not None else None).items():
                    frame.write_parquet(dest/(name+'.parquet'),compression='zstd')
            metrics.append(dict(case=case,method=method,candidate_count=ev.clonotypes.height,
                peptide_count=ev.peptides.height,alignment_count=ev.matches.height,
                **classification_metrics(scores,positives)))
        pl.DataFrame(metrics).write_csv(a.output/'metrics.tsv',separator='\t')
        case_report=dict(case=case,error_model_parameters=model_parameters,mapping_receipt_sha256=sha(mapping_path/'mapping_receipt.json') if mapping is not None else None,error_model_sha256=model_sha256,candidates=ev.clonotypes.height,peptides=ev.peptides.height,alignments=ev.matches.height,unresolved_positional_context_residues=int((~ev.background.context_valid).sum()))
        mapping_peptides=ev.peptides.select('peptide_id','seq') if a.compatibility_queries is not None else None
        del ev,observed,scores,scorer,full_scores;gc.collect()
        if a.compatibility_queries is not None:
            export_partial_mappings(a.compatibility_queries,mapping_peptides,mapping_path/'partial_compatibility',budget)
            case_report['partial_mapping_receipt_sha256']=sha(mapping_path/'partial_compatibility/mapping_receipt.json')
            del mapping_peptides;gc.collect()
        case_report['seconds']=time.perf_counter()-start
        report['cases'].append(case_report)
        (a.output/'provenance.json').write_text(json.dumps(report,indent=2)+'\n')
        print(case,time.perf_counter()-start,flush=True)

    if any(sha(repo/name)!=digest for name,digest in source.items()) or any(sha(a.input_root/name)!=digest for name,digest in raw.items()):raise ValueError('source/input drift during fresh analysis')
    if a.calibration_root is not None and any(sha(a.calibration_root/name)!=digest for name,digest in calibration_raw.items()):raise ValueError('common calibration input drift during fresh analysis')
    if a.compatibility_queries is not None and sha(a.compatibility_queries)!=partial_sha:raise ValueError('partial query drift during fresh analysis')
    report.update(status='complete',seconds=time.perf_counter()-begun,artifact_sha256={str(p.relative_to(a.output)):sha(p) for p in a.output.rglob('*') if p.is_file() and p.name!='provenance.json'})
    (a.output/'provenance.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
