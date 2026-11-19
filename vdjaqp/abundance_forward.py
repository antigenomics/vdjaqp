"""One fresh private three-acquisition AIRR draw; no inference engine or cached stages."""
import hashlib
import gzip
import shutil
import json
import math
from pathlib import Path
import subprocess
import sys
import numpy as np
import polars as pl
from .abundance_inference import REGIONS,validate_vector
from .airr_ms_inputs import prepare_airr_ms_inputs,attach_background_membership
from .airr_ms_metrics import variable_sequences,supplied_reference_regions,eligible_species,species_metrics
from .simulation_metrics import exact_peptide_placements,psm_correctness
from .readout import read_idxml
from .data_hub import sha256


def forward_draw(repertoire,target_ids,source_format,amounts,output,*,cpus=16,mssimulator=None,stage_timeout=900,expected_denominators=None,
                 min_detectability=.5,dynamic_exclusion_seconds=0.):
    """Amounts follow target CSV row order, are positive, and sum to twelve."""
    from .resources import resolve_cpu_budget
    cpus=resolve_cpu_budget(cpus)
    if not math.isfinite(min_detectability) or not 0<=min_detectability<=1:
        raise ValueError('min_detectability must be finite and in [0,1]')
    if not math.isfinite(dynamic_exclusion_seconds) or dynamic_exclusion_seconds<0:
        raise ValueError('dynamic_exclusion_seconds must be finite and nonnegative')
    candidates=pl.read_parquet(repertoire);ids=pl.read_csv(target_ids)['clone_id'].to_list()
    if candidates.height!=33052 or len(ids)!=12 or len(set(ids))!=12:raise ValueError('full search/panel population mismatch')
    candidates=attach_background_membership(candidates,inserted_reference_ids=ids)
    amounts=np.asarray(amounts,dtype=float)
    if amounts.shape!=(12,):raise ValueError('twelve explicit relative amounts required')
    output=Path(output)
    paths=prepare_airr_ms_inputs(candidates,physical_clone_ids=ids,equal_protein_amounts=False,
        protein_relative_amounts=dict(zip(ids,amounts)),source_format=source_format,output=output,
        native_protein_intensity=10000.)
    root=output.resolve();repo=Path(__file__).resolve().parents[1]
    mapping=pl.read_csv(paths['id_mapping'],separator='\t');neutral=dict(zip(mapping['clone_id'],mapping['entry_id']))
    panel=candidates.filter(pl.col('clone_id').is_in(ids)).with_columns(pl.col('clone_id').replace_strict(neutral).alias('clone_id'))
    regions=supplied_reference_regions(panel)
    if expected_denominators is not None:
        lookup={(r['chain_id'],r['region']):r['end']-r['start'] for r in regions.to_dicts() if r['region'] in REGIONS}
        actual=[[[lookup[neutral[key],region] for region in REGIONS] for key in ids]]*3
        if not np.array_equal(actual,expected_denominators):raise ValueError('reference target denominators differ before acquisition')
    physical=variable_sequences(panel).select(pl.col('clone_id').alias('chain_id'),pl.col('variable_aa').alias('protein_aa'))
    scorer_candidates=candidates.with_columns(pl.col('clone_id').replace_strict(neutral).alias('clone_id'))
    from .data import to_gene_level
    scorer_candidates=to_gene_level(scorer_candidates)
    calls=[];spectra=[];source_reports=[]
    blocks=[('solution_trypsin','Trypsin'),('gel_trypsin','Trypsin'),('shared_chymotrypsin','Chymotrypsin')]
    for label,enzyme in blocks:
        destination=root/label
        command=[sys.executable,str(repo/'scripts/smoke_proteomics.py'),str(paths['proteins_fasta'].resolve()),
            str(paths['search_fasta'].resolve()),str(destination),'--enzyme',enzyme,'--threads',str(cpus),
            '--keep-spectra','--gradient-seconds','1200','--ms2-per-bin','5','--stage-timeout',str(stage_timeout),
            '--min-detectability',str(min_detectability),'--dynamic-exclusion-seconds',str(dynamic_exclusion_seconds)]
        if mssimulator:command+=['--mssimulator',str(Path(mssimulator).resolve())]
        with (root/(label+'.log')).open('w') as log:
            subprocess.run(command,cwd=repo,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=2*stage_timeout+60)
        native=json.loads((destination/'ms2_truth.json').read_text())
        frame=read_idxml(destination/'psms.idXML',run_id=label)
        per,accuracy=psm_correctness(frame,native)
        if per.filter(pl.col('status')!='evaluated').height:raise ValueError('unlinked native assignments')
        per.write_parquet(root/(label+'_psm_correctness.parquet'),compression='zstd')
        accuracy.write_parquet(root/(label+'_accuracy.parquet'),compression='zstd')
        calls.append(frame);spectra.append(destination/'spectra.mzML');source_reports.append(command)
    queries=eligible_species(pl.concat(calls,how='vertical_relaxed')).with_row_index('index').with_columns(
        pl.col('index').cast(pl.String).alias('peptide_id')).drop('index')
    placements=exact_peptide_placements(physical,queries,threads=cpus)
    yields=[];coverages=[];denominators=[]
    for (label,_),frame in zip(blocks,calls):
        inventory,coverage,_=species_metrics(frame,queries,placements,regions)
        inventory.write_parquet(root/(label+'_species.parquet'),compression='zstd')
        coverage.write_parquet(root/(label+'_coverage.parquet'),compression='zstd')
        yields.append(inventory.filter(pl.col('quantity')=='compatible_species')['count'].item())
        lookup={(r['chain_id'],r['region']):r for r in coverage.to_dicts()}
        coverages.append([[lookup[neutral[key],region]['coverage'] for region in REGIONS] for key in ids])
        denominators.append([[lookup[neutral[key],region]['region_length'] for region in REGIONS] for key in ids])
    # Same shared chymotrypsin acquisition enters both condition unions, once per draw.
    from .pipeline import Evidence
    from .errormodel import ErrorModel,mass_prior,_E_DB_DEFAULT
    for label,tryp in [('solution',calls[0]),('gel',calls[1])]:
        peptides=pl.concat([tryp,calls[2]],how='vertical_relaxed')
        model=ErrorModel(np.array([0.,1.]),np.array([.1,.01]),mass_prior(),_E_DB_DEFAULT)
        evidence=Evidence.from_frames(scorer_candidates,peptides,error_model=model)
        scores=evidence.score();del evidence
        if scores.height!=33052 or not scores['score'].is_finite().all():raise ValueError('invalid score population')
        scores.write_parquet(root/(label+'_scores.parquet'),compression='zstd')
    summary=dict(yield_counts=yields,coverage=coverages,denominators=denominators)
    validate_vector(summary)
    (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    report=dict(status='validated_cleanup_pending',raw_input_sha256=dict(repertoire=sha256(Path(repertoire)),target_ids=sha256(Path(target_ids))),
        assay_parameters=dict(min_detectability=min_detectability,dynamic_exclusion_seconds=dynamic_exclusion_seconds),
        relative_amounts=amounts.tolist(),relative_amount_sum=float(amounts.sum()),native_total_input_intensity=120000.,
        commands=source_reports,blocks=[label for label,_ in blocks],shared_chymotrypsin_acquisitions=1,
        recovery_used_for_selection=False,source_sha256={str(p.relative_to(repo)):sha256(p) for p in
            [*repo.glob('vdjaqp/*.py'),repo/'scripts/smoke_proteomics.py',
             *[f for f in (repo/'envs/openms').rglob('*') if f.is_file()]]},
        native_reports={label:json.loads((root/label/'summary.json').read_text()) for label,_ in blocks},
        cleanup={str(p.relative_to(root)):dict(sha256=sha256(p),status='pending') for p in spectra})
    report_path=root/'joint_draw.json';report_path.write_text(json.dumps(report,indent=2)+'\n')
    for path in spectra:
        path.unlink();report['cleanup'][str(path.relative_to(root))]['status']='deleted'
        report_path.write_text(json.dumps(report,indent=2)+'\n')
    report['consumed_native_features']={}
    for label,_ in blocks:
        feature=root/label/'truth.featureXML'
        report['consumed_native_features'][str(feature.relative_to(root))]=dict(sha256=sha256(feature),bytes=feature.stat().st_size,status='pending')
        report_path.write_text(json.dumps(report,indent=2)+'\n')
        feature.unlink()
        report['consumed_native_features'][str(feature.relative_to(root))]['status']='deleted'
        report_path.write_text(json.dumps(report,indent=2)+'\n')
    report['retention']=[]
    for path in sorted(root.rglob('*')):
        if path.suffix not in ('.tsv','.idXML') or path.stat().st_size<1_000_000:continue
        raw_hash=sha256(path);packed=path.with_name(path.name+'.gz')
        with packed.open('xb') as raw:
            with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as target,path.open('rb') as source:
                shutil.copyfileobj(source,target,1024*1024)
        with gzip.open(packed,'rb') as stream:
            if hashlib.file_digest(stream,'sha256').hexdigest()!=raw_hash:raise ValueError('retention roundtrip failed')
        record=dict(raw_path=str(path.relative_to(root)),raw_sha256=raw_hash,gzip_path=str(packed.relative_to(root)),
                    gzip_sha256=sha256(packed),raw_bytes=path.stat().st_size,gzip_bytes=packed.stat().st_size,status='original_retained')
        report['retention'].append(record);report_path.write_text(json.dumps(report,indent=2)+'\n')
        path.unlink();record['status']='compressed_original_removed'
        report_path.write_text(json.dumps(report,indent=2)+'\n')
    report['status']='complete'
    report['retained_bytes']=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    report['artifact_sha256']={str(p.relative_to(root)):sha256(p) for p in root.rglob('*') if p.is_file() and p!=report_path}
    report_path.write_text(json.dumps(report,indent=2)+'\n')
    return summary
