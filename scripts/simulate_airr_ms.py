"""One fresh supplied-AIRR physical scope/acquisition case; no inference from truth.

Private outputs retain the source mapping. Each invocation rebuilds raw inputs,
simulates separate tryptic/chymotryptic acquisitions and pools evidence afterward.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time


def parse_args(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repertoire',type=Path,required=True)
    p.add_argument('--target-ids',type=Path,required=True)
    p.add_argument('--source-format',required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--physical-scope',choices=('all','panel'),required=True)
    p.add_argument('--ms2-per-bin',type=int,choices=(5,10),required=True)
    p.add_argument('--cpus',type=int)
    p.add_argument('--native-protein-intensity',type=float,default=10000.,
                   help='equal OpenMS simulator intensity units per physical record; not measured molarity')
    p.add_argument('--mssimulator',type=Path,help='explicit MSSimulator executable override')
    p.add_argument('--stage-timeout',type=float,default=900.0)
    args=p.parse_args(argv)
    if not math.isfinite(args.stage_timeout) or args.stage_timeout<=0:
        p.error('--stage-timeout must be finite and positive')
    if args.mssimulator is not None:
        args.mssimulator=args.mssimulator.resolve()
        if not args.mssimulator.is_file() or not os.access(args.mssimulator,os.X_OK):
            p.error('--mssimulator must be an executable file')
    if not math.isfinite(args.native_protein_intensity) or args.native_protein_intensity<=0:
        p.error('--native-protein-intensity must be finite and positive')
    if args.output.exists():p.error('output must be a fresh absent directory')
    if not args.repertoire.is_file() or not args.target_ids.is_file():p.error('explicit raw input files must exist')
    return args


def main(argv=None):
    args=parse_args(argv)
    repo=Path(__file__).resolve().parents[1];sys.path.insert(0,str(repo))
    from vdjaqp import resolve_cpu_budget
    cpus=resolve_cpu_budget(args.cpus)
    import numpy as np
    import polars as pl
    from vdjaqp.airr_ms_inputs import prepare_airr_ms_inputs,attach_background_membership
    from vdjaqp.airr_ms_metrics import (variable_sequences,supplied_reference_regions,eligible_species,
        species_metrics,supplied_raw_regions,variable_only_placements,supplied_gene_coverage)
    from vdjaqp.simulation_metrics import exact_peptide_placements,psm_correctness,identity_metrics
    from vdjaqp.readout import read_idxml
    from vdjaqp.pipeline import Evidence
    from vdjaqp.errormodel import ErrorModel,mass_prior,_E_DB_DEFAULT
    from vdjaqp.bench import rank_targets
    from vdjaqp.data import to_gene_level
    candidates=pl.read_parquet(args.repertoire)
    target_ids=pl.read_csv(args.target_ids)['clone_id'].to_list()
    if len(set(target_ids))!=len(target_ids) or len(target_ids)!=12:
        raise ValueError('target panel requires exactly 12 distinct supplied clone IDs')
    if not set(target_ids)<=set(candidates['clone_id']):
        raise ValueError('target panel IDs must all be present in the supplied repertoire')
    candidates=attach_background_membership(candidates,inserted_reference_ids=target_ids)
    physical_ids=candidates['clone_id'].to_list() if args.physical_scope=='all' else target_ids
    paths=prepare_airr_ms_inputs(candidates,physical_clone_ids=physical_ids,
        equal_protein_amounts=True,source_format=args.source_format,output=args.output,
        native_protein_intensity=args.native_protein_intensity)
    root=args.output.resolve()
    def save(frame,name):
        frame.write_parquet(root/(name+'.parquet'),compression='zstd')
        frame.select([name for name,dtype in frame.schema.items() if not dtype.is_nested()]).write_csv(root/(name+'.tsv'),separator='\t')
    mapping=pl.read_csv(paths['id_mapping'],separator='\t')
    # Neutral source IDs are the only identifiers entering OpenMS and scoring.
    neutral=candidates.join(mapping.select('clone_id','entry_id'),on='clone_id',validate='1:1')
    variables=variable_sequences(candidates).join(mapping.select('clone_id','entry_id'),on='clone_id',validate='1:1')
    physical=neutral.filter(pl.col('clone_id').is_in(physical_ids))
    metadata=(physical.select('entry_id',pl.col('seq').alias('protein_aa'))
        .join(variables.select('entry_id','variable_aa','variable_length'),on='entry_id',validate='1:1')
        .rename({'entry_id':'chain_id'}).with_columns(pl.lit('IGH').alias('locus')))
    panel_ids=mapping.filter(pl.col('clone_id').is_in(target_ids))['entry_id']
    panel=neutral.filter(pl.col('clone_id').is_in(target_ids))
    panel_regions=supplied_reference_regions(panel.drop('clone_id').rename({'entry_id':'clone_id'}))
    physical_regions=supplied_raw_regions(physical.drop('clone_id').rename({'entry_id':'clone_id'}))
    save(metadata,'physical_truth');save(panel_regions,'panel_regions');save(physical_regions,'physical_regions')
    scorer_candidates=to_gene_level(neutral.drop('clone_id','entry_id').with_columns(neutral['entry_id'].alias('clone_id')))
    calls=[];stages=[];spectra=[]
    for enzyme in ('Trypsin','Chymotrypsin'):
        destination=root/enzyme
        command=[sys.executable,str(repo/'scripts/smoke_proteomics.py'),str(Path(paths['proteins_fasta']).resolve()),
            str(Path(paths['search_fasta']).resolve()),str(destination),'--enzyme',enzyme,'--threads',str(cpus),
            '--keep-spectra','--gradient-seconds','1200','--ms2-per-bin',str(args.ms2_per_bin),
            '--stage-timeout',str(args.stage_timeout)]
        if args.mssimulator is not None:command.extend(['--mssimulator',str(args.mssimulator)])
        begin=time.perf_counter()
        with (root/(enzyme+'.log')).open('w') as log:
            subprocess.run(command,cwd=repo,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=3600)
        stages.append(dict(enzyme=enzyme,seconds=time.perf_counter()-begin,command=command))
        native=json.loads((destination/'ms2_truth.json').read_text())
        peptides=read_idxml(destination/'psms.idXML',run_id=enzyme).with_columns(pl.lit(enzyme).alias('enzyme'))
        per,summary=psm_correctness(peptides,native)
        if per.filter(pl.col('status')!='evaluated').height:raise ValueError('unevaluated native spectrum linkage')
        save(per,enzyme+'_psm_correctness');save(summary,enzyme+'_native_accuracy')
        save(peptides,enzyme+'_peptides');calls.append(peptides)
        spectra.append(destination/'spectra.mzML')
    pooled=pl.concat(calls,how='vertical_relaxed');save(pooled,'pooled_peptides')
    species=eligible_species(pooled)
    # One exact native batch for all eligible called species and candidate variable
    # evaluation queries. Neither source is added to the declared search FASTA.
    candidate_queries=variables.select(pl.col('variable_aa').str.replace_all('I','L',literal=True).alias('seq')).unique()
    queries=pl.concat([species,candidate_queries]).unique().sort('seq').with_row_index('index').with_columns(pl.col('index').cast(pl.String).alias('peptide_id')).drop('index')
    hits=exact_peptide_placements(metadata.select('chain_id','protein_aa'),queries,threads=cpus).join(queries,on='peptide_id',validate='m:1')
    save(hits,'exact_supplied_placements')
    variable_hits=variable_only_placements(hits,metadata)
    evaluation_queries=variables.select(pl.col('entry_id').alias('clone_id'),pl.col('variable_aa').str.replace_all('I','L',literal=True).alias('seq'))
    associations=evaluation_queries.join(variable_hits,on='seq').select('clone_id','chain_id','start','end').unique()
    reports={}
    for enzyme,peptides in [('Trypsin',calls[0]),('Chymotrypsin',calls[1]),('pooled',pooled)]:
        for scope,regions,selected_hits in [('physical',physical_regions,hits),('panel',panel_regions,variable_hits.filter(pl.col('chain_id').is_in(panel_ids.implode())))]:
            inventory,coverage,coverage_summary=species_metrics(peptides,queries,selected_hits,regions)
            for frame,suffix in [(inventory,'species'),(coverage,'coverage'),(coverage_summary,'coverage_summary')]:save(frame,enzyme+'_'+scope+'_'+suffix)
            if scope=='physical':
                save(supplied_gene_coverage(coverage),enzyme+'_gene_coverage')
        model=ErrorModel(np.array([0.,1.]),np.array([.1,.01]),mass_prior(),_E_DB_DEFAULT)
        ev=Evidence.from_frames(scorer_candidates,peptides,error_model=model)
        scores=ev.score()
        del ev
        if scores.height!=candidates.height or not np.isfinite(scores['score'].to_numpy()).all():raise ValueError('invalid score population')
        ranks=rank_targets(scores,'score',set(scores['clone_id']))
        identity,called,recovery=identity_metrics(metadata,scorer_candidates,ranks,associations)
        for frame,suffix in [(scores,'scores'),(ranks,'ranks'),(called,'candidate_identity'),(identity,'identity_population'),(recovery,'identity_metrics')]:save(frame,enzyme+'_'+suffix)
        # Stratification follows global identity resolution; never resolve an alias
        # by hiding physical competitors outside the supplied target panel.
        panel_identity=identity.filter(pl.col('chain_id').is_in(panel_ids.implode()))
        save(panel_identity,enzyme+'_panel_identity')
        classes=panel_identity.unique('variable_mass_class_id')
        panel_recovery=dict(physical_reference_rows=panel_identity.height,variable_mass_classes=classes.height,
            recovered_classes_at20=classes.filter(pl.col('class_recovered_at20')).height)
        observable=classes.filter(pl.col('full_variable_compatible_candidate_count')>0).height
        panel_recovery.update(full_variable_observable_classes=observable,
            class_recovery_at20=panel_recovery['recovered_classes_at20']/classes.height if classes.height else None,
            conditional_class_recovery_at20=panel_recovery['recovered_classes_at20']/observable if observable else None)
        save(pl.DataFrame([panel_recovery]),enzyme+'_panel_recovery')
        reports[enzyme]=dict(global_identity_metrics=recovery.to_dicts(),panel_recovery=panel_recovery)
    source_files=[repo/'scripts/simulate_airr_ms.py',repo/'scripts/smoke_proteomics.py',*repo.glob('vdjaqp/*.py'),*[p for p in (repo/'envs/openms').rglob('*') if p.is_file()]]
    def digest(path):
        with Path(path).open('rb') as stream:
            return hashlib.file_digest(stream,'sha256').hexdigest()
    report=dict(physical_scope=args.physical_scope,physical_rows=metadata.height,candidate_rows=candidates.height,
        panel_rows=12,equal_protein_amounts=True,
        native_protein_intensity_per_record=args.native_protein_intensity,
        native_protein_intensity_total=metadata.height*args.native_protein_intensity,
        coverage_conventions='physical full supplied/variable/actual Cstub and contiguous raw label runs; panel corrected FR/CDR regions for fit; label6 does not infer D/N boundaries; D is templated',
        intensity_units='OpenMS MSSimulator native protein intensity units; equal per record; not measured molarity',
        intensity_default_source='OpenMS MSSimulator.cpp lines203-230 commit60ce8966d577b43a6419403dec171539e08fe513',
        gradient_seconds=1200,stage_timeout_seconds=args.stage_timeout,
        mssimulator_override=str(args.mssimulator) if args.mssimulator else None,ms2_per_bin=args.ms2_per_bin,
        cpus=cpus,stages=stages,readouts=reports,threshold_confirmed_nats=20,
        raw_input_sha256=dict(repertoire=digest(args.repertoire),target_ids=digest(args.target_ids)),
        software_sha256={str(p.relative_to(repo)):digest(p) for p in source_files},
        artifact_sha256={str(p.relative_to(root)):digest(p) for p in root.rglob('*') if p.is_file() and p not in spectra},
        spectra_cleanup={str(p.relative_to(root)):dict(sha256=digest(p),bytes=p.stat().st_size,status='pending') for p in spectra})
    # All validation and retained tables precede removal. Failed runs keep spectra.
    (root/'airr_ms_report.json').write_text(json.dumps(report,indent=2)+'\n')
    for path in spectra:
        path.unlink()
        report['spectra_cleanup'][str(path.relative_to(root))]['status']='deleted'
        (root/'airr_ms_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(reports,indent=2))


if __name__=='__main__':
    main()
