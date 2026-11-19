"""Fresh role-correct scoring of eighteen retained AIRR contexts; no native MS.

Read raw call fields from final PSM exports, never fitted alignment/features or
truth annotations. Accepted historical scores remain immutable comparators.
"""
import argparse,gc,hashlib,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

CASE_NAMES={'baseline','exclusion30','detect075','detect075_exclusion30','detect085','detect090'}
SETTINGS={'baseline':(.5,0.),'exclusion30':(.5,30.),'detect075':(.75,0.),
          'detect075_exclusion30':(.75,30.),'detect085':(.85,0.),'detect090':(.9,0.)}
BLOCKS=('solution_trypsin','gel_trypsin','shared_chymotrypsin')
STUDY_NAMES={'panel5','panel10','all5','all10','loading2500','loading5000'}
STUDY_SETTINGS={'panel5':('panel',5,10000.),'panel10':('panel',10,10000.),
    'all5':('all',5,10000.),'all10':('all',10,10000.),
    'loading2500':('panel',5,2500.),'loading5000':('panel',5,5000.)}


def study_calls(frame):
    """Keep both original protease identities, stripping any non-readout fields."""
    import polars as pl
    if 'enzyme' not in frame.columns or set(frame['enzyme'])!={'Trypsin','Chymotrypsin'}:
        raise ValueError('original pooled calls require both named proteases')
    parts=[]
    for enzyme in ['Trypsin','Chymotrypsin']:
        parts.append(raw_calls(frame.filter(pl.col('enzyme')==enzyme),enzyme).with_columns(pl.lit(enzyme).alias('enzyme')))
    return pl.concat(parts,how='vertical_relaxed')


def scoped_metrics(frame,physical,references):
    """Separate physical presence from reference retrieval in all-physical runs."""
    from vdjaqp.bench import classification_metrics
    row=classification_metrics(frame,physical)
    row.update(metric_scope='physical_record_presence',positive_definition='physically simulated candidate records',
               negative_definition='candidate records absent from declared physical input')
    if len(physical)==frame.height:
        row.update(average_precision=None,roc_auc=None,negative_definition='none: all candidate records physically present; AP/AUROC unavailable')
    result=[row]
    if len(physical)==frame.height:
        retrieval=classification_metrics(frame,references)
        retrieval.update(metric_scope='reference_retrieval',positive_definition='12 inserted reference records',
            negative_definition='nonreference candidate records; physically present, not identification negatives')
        result.append(retrieval)
    return result


def raw_calls(frame,run_id):
    """Whitelist published readout fields, stripping every evaluation/truth column."""
    import polars as pl
    from vdjaqp.readout import _SCHEMA
    if not set(_SCHEMA)<=set(frame.columns):raise ValueError('missing raw readout fields')
    calls=frame.select(list(_SCHEMA)).cast(_SCHEMA,strict=True)
    if (calls['peptide_id'].n_unique()!=calls.height or calls['peptide_id'].null_count() or calls['spectrum_id'].null_count()
            or calls.filter(~pl.col('peptide_id').str.starts_with(run_id+':psm:') |
                            ~pl.col('spectrum_id').str.starts_with(run_id+':')).height):
        raise ValueError('raw call IDs do not match original acquisition run_id')
    return calls


def mapped_candidates(candidates,mapping,references):
    """Apply the recorded bijection after attaching source, not physical, role."""
    import polars as pl
    from vdjaqp.airr_ms_inputs import attach_background_membership
    if (set(mapping['clone_id'])!=set(candidates['clone_id']) or
            mapping['clone_id'].n_unique()!=mapping.height or mapping['entry_id'].n_unique()!=mapping.height or
            mapping['clone_id'].null_count() or mapping['entry_id'].null_count()):
        raise ValueError('mapping must be a complete candidate bijection')
    tagged=attach_background_membership(candidates,inserted_reference_ids=references)
    if 'background_member' in mapping.columns:
        if mapping.schema['background_member']!=pl.Boolean or mapping['background_member'].null_count():
            raise ValueError('mapping background_member must be Boolean and nonmissing')
        check=mapping.select('clone_id','background_member').join(tagged.select('clone_id','background_member'),
            on='clone_id',suffix='_source',validate='1:1')
        if check.filter(pl.col('background_member')!=pl.col('background_member_source')).height:
            raise ValueError('mapping background role conflicts with source')
    lookup=dict(zip(mapping['clone_id'],mapping['entry_id']))
    return tagged.with_columns(pl.col('clone_id').replace_strict(lookup).alias('clone_id'))


def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root',type=Path,required=True);p.add_argument('--cases',type=Path,required=True,
        help='JSON list of six {case,path} records pointing to completed joint draws')
    p.add_argument('--output',type=Path,required=True);p.add_argument('--cpus',type=int,default=16)
    p.add_argument('--study-cases',type=Path,required=True,help='JSON list of six original {case,path} study records')
    a=p.parse_args()
    from vdjaqp import resolve_cpu_budget
    budget=resolve_cpu_budget(a.cpus)
    specs=json.loads(a.cases.read_text())
    study_specs=json.loads(a.study_cases.read_text())
    if len(specs)!=6 or {x['case'] for x in specs}!=CASE_NAMES:p.error('exact six-condition inventory required')
    if len(study_specs)!=6 or {x['case'] for x in study_specs}!=STUDY_NAMES:p.error('exact six original study cases required')
    if a.output.exists():p.error('fresh output required')
    import numpy as np
    import polars as pl
    from vdjaqp.pipeline import Evidence
    from vdjaqp.errormodel import ErrorModel,mass_prior
    from vdjaqp.bench import classification_metrics,rank_targets
    from vdjaqp.data import to_gene_level
    repo=Path(__file__).resolve().parents[1]
    repertoire=a.input_root/'dataset_mixa_gel/clonotypes.parquet';target=a.input_root/'ground_truth/target_clonotypes.csv'
    candidates=pl.read_parquet(repertoire);references=pl.read_csv(target)['clone_id'].to_list()
    if candidates.height!=33052 or candidates['clone_id'].n_unique()!=33052 or len(set(references))!=12 or len(references)!=12:raise ValueError('source population mismatch')
    sources={str(f.relative_to(repo)):sha(f) for f in [*sorted((repo/'vdjaqp').glob('*.py')),Path(__file__).resolve()]}
    study_repertoire=a.input_root/'dataset_mixa/clonotypes.parquet'
    inputs={str(f.resolve()):sha(f) for f in [repertoire,study_repertoire,target,a.cases,a.study_cases]}
    a.output.mkdir(parents=True);metrics=[];deltas=[]
    report=dict(status='running',source_sha256=sources,input_sha256=inputs,cpu_budget=budget,
        score_model='fixed ErrorModel endpoints .1/.01; n_substitutions=0; e_db=.1',
        aggregation='alignment_union',new_native_acquisitions=0,recovery_used_for_assay_selection=False,cases=[])
    def save(): (a.output/'provenance.json').write_text(json.dumps(report,indent=2)+'\n')
    save()
    try:
        for spec in specs:
            case=spec['case'];root=Path(spec['path']).resolve();joint=json.loads((root/'joint_draw.json').read_text())
            if joint['status']!='complete' or any(joint['raw_input_sha256'][k]!=sha(f) for k,f in [('repertoire',repertoire),('target_ids',target)]):
                raise ValueError('case source input/status differs')
            settings=joint.get('assay_parameters',dict(min_detectability=.5,dynamic_exclusion_seconds=0.))
            if (tuple(settings[k] for k in ['min_detectability','dynamic_exclusion_seconds'])!=SETTINGS[case]
                    or joint['relative_amounts']!=[1.]*12):
                raise ValueError('case label/settings or equal-source amounts differ')
            mapping=pl.read_csv(root/'private_id_mapping.tsv',separator='\t')
            clones=to_gene_level(mapped_candidates(candidates,mapping,references))
            if clones['background_member'].sum()!=33040:raise ValueError('require 33040 background /12 inserted source roles')
            if mapping.schema['physical']!=pl.Boolean or mapping['physical'].null_count():raise ValueError('invalid physical evaluation flag')
            physical=set(mapping.filter('physical')['entry_id'])
            if len(physical)!=12:raise ValueError('joint assay physical evaluation population must be12')
            paths=[root/'joint_draw.json',root/'private_id_mapping.tsv',root/'summary.json',root/'discrepancy.json',
                   *[root/(b+'_psm_correctness.parquet') for b in BLOCKS],*[root/(x+'_scores.parquet') for x in ['solution','gel']]]
            hashes={str(f):sha(f) for f in paths}
            for f in paths:
                relative=str(f.relative_to(root))
                if relative in joint.get('artifact_sha256',{}) and joint['artifact_sha256'][relative]!=hashes[str(f)]:
                    raise ValueError('retained artifact hash differs')
            calls=[raw_calls(pl.read_parquet(root/(b+'_psm_correctness.parquet')),b) for b in BLOCKS]
            dest=a.output/case;dest.mkdir()
            for condition,index in [('solution',0),('gel',1)]:
                observed=pl.concat([calls[index],calls[2]],how='vertical_relaxed')
                model=ErrorModel(np.array([0.,1.]),np.array([.1,.01]),mass_prior(),0,.1)
                evidence=Evidence.from_frames(clones,observed,error_model=model)
                scores=evidence.score();accepted=pl.read_parquet(root/(condition+'_scores.parquet'))
                for frame in [scores,accepted]:
                    if frame.height!=33052 or set(frame['clone_id'])!=set(clones['clone_id']):raise ValueError('candidate population changed')
                scores.write_parquet(dest/(condition+'_scores.parquet'),compression='zstd')
                rank_targets(scores,'score',physical).write_parquet(dest/(condition+'_target_ranks.parquet'),compression='zstd')
                rank_targets(accepted,'score',physical).write_parquet(dest/(condition+'_accepted_target_ranks.parquet'),compression='zstd')
                for method,frame in [('accepted_erased_role',accepted),('proper_source_role',scores)]:
                    metrics.extend(dict(case=case,condition=condition,method=method,**row)
                                   for row in scoped_metrics(frame,physical,physical))
                x=accepted.sort('clone_id');y=scores.sort('clone_id')
                for column in ['confirmed_nats','junction_nats','contradiction','evidence','score','framework_only','n_tied']:
                    d=y[column].cast(pl.Float64).to_numpy()-x[column].cast(pl.Float64).to_numpy()
                    deltas.append(dict(case=case,condition=condition,quantity=column,n_candidates=33052,
                        n_exactly_changed=int(np.count_nonzero(d)),min_delta=float(d.min()),max_delta=float(d.max())))
                pl.DataFrame(metrics).write_csv(a.output/'metrics.tsv',separator='\t')
                pl.DataFrame(deltas).write_csv(a.output/'candidate_deltas.tsv',separator='\t')
                print(json.dumps(dict(case=case,condition=condition,candidates=scores.height,peptides=observed.height,alignments=evidence.alignments.height)),flush=True)
                del evidence,scores,accepted;gc.collect()
            if any(sha(Path(f))!=h for f,h in hashes.items()):raise ValueError('case input artifact drift')
            report['cases'].append(dict(case=case,input_sha256=hashes,background_members=33040,inserted_references=12,
                physical_evaluation_records=12,assay_parameters=settings,summary_and_discrepancy_unchanged=True));save()
        study_candidates=pl.read_parquet(study_repertoire)
        from vdjaqp.airr_ms_metrics import variable_sequences,variable_only_placements
        from vdjaqp.simulation_metrics import identity_metrics
        for spec in study_specs:
            case='study_'+spec['case'];root=Path(spec['path']).resolve()
            producer=root/'airr_ms_report.json';prior=json.loads(producer.read_text())
            if (prior['physical_scope'],prior['ms2_per_bin'],prior['native_protein_intensity_per_record'])!=STUDY_SETTINGS[spec['case']]:
                raise ValueError('original study case/settings differ')
            if prior['raw_input_sha256']!={'repertoire':sha(study_repertoire),'target_ids':sha(target)}:
                raise ValueError('original study raw source differs')
            mapping_path=root/'private_id_mapping.tsv'
            if not mapping_path.exists():mapping_path=root/'private_id_mapping.tsv.gz'
            mapping=pl.read_csv(mapping_path,separator='\t');clones=to_gene_level(mapped_candidates(study_candidates,mapping,references))
            if clones.height!=33052 or clones['background_member'].sum()!=33040:raise ValueError('study source-role population differs')
            physical=set(mapping.filter('physical')['entry_id']);ref_entries=set(mapping.filter(pl.col('clone_id').is_in(references))['entry_id'])
            if len(physical)!=prior['physical_rows'] or len(ref_entries)!=12:raise ValueError('study physical/reference scope differs')
            paths=[producer,mapping_path,*[root/n for n in ['pooled_peptides.parquet','pooled_scores.parquet',
                'physical_truth.parquet','exact_supplied_placements.parquet','pooled_identity_metrics.parquet']]]
            hashes={str(f):sha(f) for f in paths}
            for f in paths[1:]:
                name=f.name;expected=prior['artifact_sha256'].get(name)
                if name=='private_id_mapping.tsv.gz':
                    import gzip
                    with gzip.open(f,'rb') as stream:actual=hashlib.file_digest(stream,'sha256').hexdigest()
                    expected=prior['artifact_sha256'].get('private_id_mapping.tsv')
                else:actual=hashes[str(f)]
                if expected is None or actual!=expected:raise ValueError('original study retained artifact hash mismatch')
            observed=study_calls(pl.read_parquet(root/'pooled_peptides.parquet'))
            model=ErrorModel(np.array([0.,1.]),np.array([.1,.01]),mass_prior(),0,.1)
            evidence=Evidence.from_frames(clones,observed,error_model=model);scores=evidence.score()
            accepted=pl.read_parquet(root/'pooled_scores.parquet')
            for frame in [scores,accepted]:
                if frame.height!=33052 or set(frame['clone_id'])!=set(clones['clone_id']):raise ValueError('study candidate population differs')
            dest=a.output/case;dest.mkdir();scores.write_parquet(dest/'pooled_scores.parquet',compression='zstd')
            # Final exact placements are evaluation-only; never enter Evidence.
            metadata=pl.read_parquet(root/'physical_truth.parquet')
            if set(metadata['chain_id'])!=physical:raise ValueError('physical truth and explicit membership differ')
            placements=variable_only_placements(pl.read_parquet(root/'exact_supplied_placements.parquet'),metadata)
            variables=variable_sequences(clones).select('clone_id',pl.col('variable_aa').str.replace_all('I','L',literal=True).alias('seq'))
            associations=variables.join(placements,on='seq').select('clone_id','chain_id','start','end').unique()
            summaries=[]
            for method,frame in [('accepted_erased_role',accepted),('proper_source_role',scores)]:
                metrics.extend(dict(case=case,condition='pooled',method=method,**row) for row in scoped_metrics(frame,physical,ref_entries))
                rank_targets(frame,'score',ref_entries).write_parquet(dest/(method+'_reference_ranks.parquet'),compression='zstd')
                ranks=rank_targets(frame,'score',set(frame['clone_id']))
                population,called,summary=identity_metrics(metadata,clones,ranks,associations)
                summaries.append(summary.with_columns(pl.lit(method).alias('method')))
                population.write_parquet(dest/(method+'_identity_population.parquet'),compression='zstd')
                called.write_parquet(dest/(method+'_candidate_identity.parquet'),compression='zstd')
                if method=='accepted_erased_role':
                    saved=pl.read_parquet(root/'pooled_identity_metrics.parquet')
                    if not summary.equals(saved):raise ValueError('accepted identity metrics fail exact reconstruction')
            pl.concat(summaries).write_csv(dest/'identity_metrics.tsv',separator='\t')
            x=accepted.sort('clone_id');y=scores.sort('clone_id')
            for column in ['confirmed_nats','junction_nats','contradiction','evidence','score','framework_only','n_tied']:
                d=y[column].cast(pl.Float64).to_numpy()-x[column].cast(pl.Float64).to_numpy()
                deltas.append(dict(case=case,condition='pooled',quantity=column,n_candidates=33052,
                    n_exactly_changed=int(np.count_nonzero(d)),min_delta=float(d.min()),max_delta=float(d.max())))
            pl.DataFrame(metrics).write_csv(a.output/'metrics.tsv',separator='\t');pl.DataFrame(deltas).write_csv(a.output/'candidate_deltas.tsv',separator='\t')
            if any(sha(Path(f))!=h for f,h in hashes.items()):raise ValueError('study input artifact drift')
            report['cases'].append(dict(case=case,condition='pooled',input_sha256=hashes,
                background_members=33040,inserted_references=12,physical_evaluation_records=len(physical),
                exact_final_placements_used_only_for_evaluation=True));save()
            print(json.dumps(dict(case=case,condition='pooled',candidates=scores.height,alignments=evidence.alignments.height)),flush=True)
            del evidence,scores,accepted;gc.collect()
        if any(sha(repo/f)!=h for f,h in sources.items()) or any(sha(Path(f))!=h for f,h in inputs.items()):raise ValueError('source/raw drift')
        report['status']='complete'
    except BaseException as error:report.update(status='failed',error=repr(error));raise
    finally:
        report['artifact_sha256']={str(f.relative_to(a.output)):sha(f) for f in a.output.rglob('*') if f.is_file() and f.name!='provenance.json'};save()


if __name__=='__main__':main()
