"""Final observation-only evaluation; truth never constructs inference inputs.

One all-chain annotation batch and one deduplicated observed-AA placement batch.
No MS, peptide scoring, fake ranks or identification recovery are performed.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def raw_truth_digest(path):
    h=hashlib.sha256();n=0
    with gzip.open(path,'rb') as stream:
        while data:=stream.read(1024*1024):h.update(data);n+=len(data)
    return h.hexdigest(),n


def queries_and_aliases(entries,candidates):
    import polars as pl
    provenance=['sequence_origin','coverage_kind','source_observation_id','observed_read_id','parent_clone_id','native_clone_row','duplicate_count','native_junction_nt','native_junction_aa']
    frames=[]
    for frame,key,level in [(entries,'entry_id','arda_observed_exact_compatibility'),(candidates,'clone_id','heavy_candidate')]:
        frames.append(frame.select(pl.col(key).alias('observed_id'),'seq',*[c for c in provenance if c in frame.columns]).with_columns(pl.lit(level).alias('evidence_level')))
    aliases=pl.concat(frames,how='diagonal_relaxed')
    valid=pl.col('seq').is_not_null() & pl.col('seq').str.contains('^[ACDEFGHIKLMNPQRSTVWY]{6,}$')
    dropped=aliases.filter(~valid.fill_null(False))
    queries=aliases.filter(valid).select('seq').unique().sort('seq').with_row_index('_id').with_columns(pl.format('observed_{}',pl.col('_id')).alias('peptide_id')).drop('_id')
    return queries,aliases.filter(valid).join(queries,on='seq',how='left',validate='m:1'),dropped


def observability(metadata,hits,aliases):
    """Same pair/full criterion as identity_metrics, without invented scores."""
    import polars as pl
    from vdjaqp.match import ms_collapse
    heavy=metadata.filter(pl.col('locus')=='IGH')
    labels=['mass_'+hashlib.sha256(seq.encode()).hexdigest() for seq in ms_collapse(heavy['variable_aa'])]
    heavy=heavy.with_columns(pl.Series('variable_mass_class_id',labels,dtype=pl.String))
    pairs=(hits.join(aliases.filter(pl.col('evidence_level')=='heavy_candidate').select('peptide_id',pl.col('observed_id').alias('clone_id')),on='peptide_id')
           .join(heavy.select('chain_id','variable_mass_class_id','variable_length'),on='chain_id')
           .with_columns(((pl.col('start')==0)&(pl.col('end')>=pl.col('variable_length'))).alias('full_variable_compatible')))
    pairs=pairs.with_columns((pl.min_horizontal('end','variable_length')-pl.col('start')).clip(lower_bound=0).alias('variable_overlap_residues'))
    associations=pairs.select('clone_id','chain_id','start','end','full_variable_compatible','variable_overlap_residues').unique()
    classes=(pairs.group_by('clone_id','variable_mass_class_id').agg(pl.col('full_variable_compatible').any(),(pl.col('variable_overlap_residues')>0).any().alias('_variable_overlap')))
    totals=classes.group_by('variable_mass_class_id').agg(pl.col('full_variable_compatible').sum().alias('full_variable_compatible_candidate_count'),((~pl.col('full_variable_compatible'))&pl.col('_variable_overlap')).sum().alias('partial_compatible_candidate_count'))
    return heavy.join(totals,on='variable_mass_class_id',how='left').with_columns(pl.col('full_variable_compatible_candidate_count','partial_compatible_candidate_count').fill_null(0)),associations



def associated_gene_capture(metadata,supported_hits):
    """Any supported protein placement, distinct from direct segment coverage."""
    import polars as pl
    supported=supported_hits.select('chain_id').unique().with_columns(pl.lit(True).alias('supported'))
    m=metadata.join(supported,on='chain_id',how='left').with_columns(pl.col('supported').fill_null(False),pl.col('variable_aa').str.replace_all('I','L').alias('_variable_class'))
    frames=[]
    for seg,column in [('V','v_call'),('J','j_call'),('C','constant_nt_gene')]:
        f=m.with_columns(pl.col(column).fill_null('missing').alias('gene'),pl.lit('native_allele' if seg!='C' else 'prescribed_physical_constant_gene').alias('namespace'))
        if seg=='V':
            unique=(pl.col('locus')=='IGH')&(pl.col('canonical_v_match_status')=='unique_gene')
            f=f.with_columns(pl.when(unique).then(pl.col('canonical_v_gene')).otherwise(pl.col('gene')).alias('gene'),pl.when(unique).then(pl.lit('exact_mapped_canonical_gene')).otherwise(pl.col('namespace')).alias('namespace'))
        frames.append(f.group_by('locus','gene','namespace').agg(pl.len().alias('physical_records'),pl.col('supported').sum().alias('supported_records'),pl.col('_variable_class').n_unique().alias('physical_variable_IL_classes'),pl.col('_variable_class').filter(pl.col('supported')).n_unique().alias('supported_variable_IL_classes')).with_columns(pl.lit(seg).alias('segment')))
    return pl.concat(frames).sort('locus','segment','gene')


def bounded_coverage(regions,placements,levels,*,chains_per_chunk=1024,max_slots=2000000):
    """Whole-chain blocks bound existing native table residue expansion."""
    import polars as pl
    from vdjaqp.simulation_metrics import region_coverage
    if chains_per_chunk<1 or max_slots<1:raise ValueError('positive coverage bounds required')
    blocks=regions.select('chain_id').unique().sort('chain_id').with_row_index('_index').with_columns((pl.col('_index')//chains_per_chunk).alias('_chunk')).drop('_index')
    rg=regions.join(blocks,on='chain_id');ps=placements.join(blocks,on='chain_id',how='left')
    if ps['_chunk'].null_count():raise ValueError('unknown placement chain')
    cost=pl.concat([rg.select('_chunk',((pl.col('end')-pl.col('start')).fill_null(0)*len(levels)).alias('_slots')),ps.select('_chunk',(pl.col('end')-pl.col('start')).cast(pl.Int64).alias('_slots'))]).group_by('_chunk').agg(pl.col('_slots').sum())
    if cost.filter(pl.col('_slots')>max_slots).height:raise ValueError('coverage expansion exceeds whole-chain chunk bound')
    parts=[]
    for chunk in blocks['_chunk'].unique().sort():
        parts.append(region_coverage(rg.filter(pl.col('_chunk')==chunk).drop('_chunk'),ps.filter(pl.col('_chunk')==chunk).drop('_chunk'),levels))
    return pl.concat(parts).sort('chain_id','region','evidence_level'),cost


def verify_observations(receipt_path,geometry):
    """Verify canonical native ARDA export artifacts and source provenance."""
    record=json.loads(receipt_path.read_text());base=receipt_path.parent.resolve()
    if record.get('status')!='complete' or record.get('source_format')!='arda_fastq_native' or record.get('geometry')!=geometry:
        raise ValueError('incomplete/wrong native ARDA observation receipt')
    artifacts=record['artifact_sha256']
    for name,h in artifacts.items():
        path=(base/name).resolve()
        if not path.is_relative_to(base) or digest(path)!=h:raise ValueError('observed artifact hash mismatch')
    if not {'catalog.parquet','candidates.parquet','search_entries.parquet'}<=set(artifacts):raise ValueError('missing ARDA final observation artifact')
    import polars as pl
    for name,key in [('catalog','source_observation_id'),('candidates','clone_id'),('search_entries','entry_id')]:
        frame=pl.read_parquet(base/(name+'.parquet'),columns=[key,'background_member'])
        if frame.height!=record['counts'][name] or frame[key].null_count() or frame[key].n_unique()!=frame.height or frame['background_member'].dtype!=pl.Boolean or frame['background_member'].null_count() or not frame['background_member'].all():raise ValueError('native ARDA count/key/background-role mismatch')
    checked=[];unavailable={}
    for name,h in record.get('source_input_sha256',{}).items():
        path=Path(name)
        # Collected remote source pins remain provenance, not falsely missing local inputs.
        if path.is_file():
            if digest(path)!=h:raise ValueError('native ARDA source hash mismatch')
            checked.append(name)
        else:unavailable[name]=h
    result=dict(record,source_verification=dict(accessible_source_paths_verified=checked,unavailable_source_pins=unavailable,scope='export artifacts verified locally; unavailable remote source pins retained without a byte-verification claim'))
    return result,result,base


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['truth-gzip','validation-receipt','observed-receipt','output']:parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--geometry',choices=['amplicon','rnaseq'],required=True)
    parser.add_argument('--cpus',type=int,default=16)
    parser.add_argument('--max-placements',type=int,default=5000000)
    parser.add_argument('--max-expanded-slots',type=int,default=2000000)
    args=parser.parse_args(argv)
    if args.output.exists():raise FileExistsError('fresh output required')
    if args.max_placements<1 or args.max_expanded_slots<1:raise ValueError('positive memory bounds required')
    from vdjaqp import resolve_cpu_budget
    cpus=resolve_cpu_budget(args.cpus)
    import polars as pl
    from vdjaqp.simulation_truth import annotate_truth
    from vdjaqp.simulation_metrics import exact_peptide_placements
    from scripts.smoke_metrics import aggregate_coverage,aggregate_gene_coverage
    started=time.monotonic();pins={str(p.resolve()):digest(p) for p in [args.truth_gzip,args.validation_receipt,args.observed_receipt]}
    validation=json.loads(args.validation_receipt.read_text())
    if validation['status']!='complete':raise ValueError('incomplete native validation')
    raw,rawbytes=raw_truth_digest(args.truth_gzip)
    if (raw,rawbytes)!=(validation['truth']['raw_sha256'],validation['truth']['raw_bytes']) or digest(args.truth_gzip)!=validation['truth']['gzip_sha256']:raise ValueError('truth receipt hash mismatch')
    _,receipt,observed=verify_observations(args.observed_receipt,args.geometry)
    tables={n:pl.read_parquet(observed/(n+'.parquet')) for n in ['catalog','candidates','search_entries']}
    for n,f in tables.items():
        if f.height!=receipt['counts'][n] or f['background_member'].null_count() or not f['background_member'].all():raise ValueError('observed membership/count mismatch')
    for name,key in [('catalog','source_observation_id'),('candidates','clone_id'),('search_entries','entry_id')]:
        if tables[name][key].null_count() or tables[name][key].n_unique()!=tables[name].height:raise ValueError('observed keys are not unique')
    with gzip.open(args.truth_gzip,'rt') as stream:truth=json.load(stream)
    chains=truth['chains']
    if len(chains)!=validation['observed_chains'] or len({r['chain_id'] for r in chains})!=len(chains):raise ValueError('truth population mismatch')
    if truth['heavy_recombination_model']['sha256']!=validation['model_sha256']:raise ValueError('native model identity differs')
    from collections import Counter
    if dict(Counter(c['locus'] for c in chains))!=validation['locus_counts']:raise ValueError('native locus population differs')
    q,aliases,dropped=queries_and_aliases(tables['search_entries'],tables['candidates'])
    # Exactly one native all-chain annotation, then one deduplicated AA search.
    annotation_started=time.monotonic()
    metadata,regions=annotate_truth(chains,threads=cpus);del truth,chains
    annotation_seconds=time.monotonic()-annotation_started
    before=time.monotonic();hits=exact_peptide_placements(metadata.select('chain_id','protein_aa'),q,threads=cpus);search_seconds=time.monotonic()-before
    if hits.height>args.max_placements:raise ValueError('native placement count exceeds explicit memory bound')
    level='arda_observed_exact_compatibility'
    eligible=aliases.filter(pl.col('evidence_level')==level).select('peptide_id').unique()
    placements=hits.join(eligible,on='peptide_id',how='semi').select('chain_id','start','end').unique().with_columns(pl.lit(level).alias('evidence_level'))
    coverage,cost=bounded_coverage(regions,placements,[level],max_slots=args.max_expanded_slots)
    candidate_counts=aliases.filter(pl.col('evidence_level')=='heavy_candidate').group_by('peptide_id').len()
    association_estimate=hits.join(candidate_counts,on='peptide_id').select(pl.col('len').sum()).item() or 0
    if association_estimate>args.max_placements:raise ValueError('candidate alias association expansion exceeds explicit bound')
    population,associations=observability(metadata,hits,aliases)
    classes=population.unique('variable_mass_class_id')
    summary=dict(true_heavy_rows=population.height,true_variable_mass_classes=classes.height,full_variable_observable_classes=classes.filter(pl.col('full_variable_compatible_candidate_count')>0).height,partial_only_classes=classes.filter((pl.col('full_variable_compatible_candidate_count')==0)&(pl.col('partial_compatible_candidate_count')>0)).height)
    args.output.mkdir(parents=True)
    outputs={'truth_metadata':metadata,'truth_regions':regions,'queries':q,'observed_aliases':aliases,'excluded_queries':dropped,'observed_exact_placements':hits,'chain_region_coverage':coverage,'gene_coverage':aggregate_gene_coverage(coverage,metadata),'gene_associated_capture':associated_gene_capture(metadata,hits.join(eligible,on='peptide_id',how='semi')),'locus_coverage':aggregate_coverage(coverage,['locus','region','evidence_level']),'heavy_identity_population':population,'candidate_truth_associations':associations,'coverage_expansion_profile':cost,'catalog_status_counts':tables['catalog'].group_by('locus','status').len()}
    for n,f in outputs.items():f.write_parquet(args.output/(n+'.parquet'))
    verify_observations(args.observed_receipt,args.geometry)
    if any(digest(Path(p))!=h for p,h in pins.items()):raise ValueError('inputs changed during evaluation')
    report=dict(status='complete',geometry=args.geometry,cpu_budget=cpus,source_input_sha256=pins,observed_counts=receipt['counts'],observability=summary,unique_queries=q.height,query_aliases=aliases.height,placement_rows=hits.height,native_search_seconds=search_seconds,native_annotation_seconds=annotation_seconds,coverage_expanded_slots=int(cost['_slots'].sum()),max_placements=args.max_placements,max_expanded_slots=args.max_expanded_slots,candidate_alias_association_upper_bound=association_estimate,observed_native_receipt=receipt,scope='evaluation-only exact compatibility; no MS/recovery/accepted identification; no truth enters observed frames',artifact_sha256={f.name:digest(f) for f in args.output.iterdir()},wall_seconds=time.monotonic()-started)
    report['source_sha256']={str(Path(sys.modules[n].__file__).resolve()):digest(Path(sys.modules[n].__file__)) for n in [__name__,'vdjaqp.simulation_truth','vdjaqp.simulation_metrics','scripts.smoke_metrics']}
    (args.output/'provenance.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
