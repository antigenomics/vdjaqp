"""All-locus observed compatibility capture; no light scoring or identity inference.

One exact I/L-compatible batch over distinct exported native target sequences and
all called peptide AA. Existing preparation pools share one chymotrypsin input.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.export_mock_observation_metrics import bounded_coverage,digest
from scripts.run_matched_reference_ms import BLOCKS,verify_manifest

CONDITIONS=('solution','heavy_band')
UNIT='native_processed_fragment_support'
LEVEL='called_peptide_exact_IL_compatibility'


def load_observed(receipt_path,native_clones,geometry):
    """Verify final frames and the supplied exact native clone source, not other raw TSVs."""
    import polars as pl
    receipt_path,native_clones=Path(receipt_path),Path(native_clones)
    receipt=json.loads(receipt_path.read_text());base=receipt_path.parent.resolve()
    if (receipt.get('status'),receipt.get('source_format'),receipt.get('geometry'))!=('complete','arda_fastq_native',geometry):raise ValueError('complete matching native ARDA receipt required')
    frames={}
    for name,key in [('catalog','source_observation_id'),('candidates','clone_id'),('search_entries','entry_id')]:
        p=base/(name+'.parquet')
        if digest(p)!=receipt['artifact_sha256'][p.name]:raise ValueError('observed artifact hash mismatch')
        f=pl.read_parquet(p)
        if not {key,'background_member'}<=set(f.columns):raise ValueError('required observed schema fields missing')
        if f.height!=receipt['counts'][name] or f[key].null_count() or f[key].n_unique()!=f.height or f['background_member'].dtype!=pl.Boolean or f['background_member'].null_count() or not f['background_member'].all():raise ValueError('observed key/count/background-role mismatch')
        frames[name]=f
    expected={h for p,h in receipt['source_input_sha256'].items() if Path(p).name.endswith('.clones.tsv')}
    if len(expected)!=1 or digest(native_clones) not in expected:raise ValueError('native clone TSV must match producing source hash')
    catalog,entries=frames['catalog'],frames['search_entries']
    required={'source_observation_id','native_clone_row','parent_clone_id','locus','native_junction_aa','duplicate_count'}
    if not (required|{'coverage_kind','sequence_origin'})<=set(catalog.columns) or not {'source_observation_id','seq'}<=set(entries.columns) or (entries.height and not required<=set(entries.columns)):raise ValueError('required native observation schema fields missing')
    if not catalog.schema['native_clone_row'].is_integer() or not catalog.schema['duplicate_count'].is_integer() or catalog['duplicate_count'].null_count() or catalog['native_clone_row'].null_count() or catalog['native_clone_row'].n_unique()!=catalog.height or catalog.filter(~pl.col('locus').is_in(['IGH','IGK','IGL'])).height:raise ValueError('invalid all-locus native catalog population')
    raw=pl.scan_csv(native_clones,separator='\t',quote_char=None,infer_schema_length=0,row_index_name='native_clone_row')
    columns=['locus','junction_aa','v_call','j_call','c_call','duplicate_count']
    if not set(columns)<=set(raw.collect_schema().names()):raise ValueError('native clone metadata columns missing')
    clones=raw.select(pl.col('native_clone_row').cast(pl.Int64),*[pl.when(pl.col(c)=='""').then(None).otherwise(pl.col(c)).alias(c) for c in columns]).collect(engine='streaming')
    clones=clones.with_columns(pl.col('duplicate_count').cast(pl.Int64,strict=True))
    if clones.filter(pl.col('duplicate_count').is_null()|(pl.col('duplicate_count')<0)).height:raise ValueError('native support counts must be nonnegative integers')
    records=catalog.join(clones,on='native_clone_row',how='left',suffix='_source',validate='1:1')
    if records.filter(pl.col('locus_source').is_null()|(pl.col('locus')!=pl.col('locus_source'))|(pl.col('native_junction_aa')!=pl.col('junction_aa'))|(pl.col('duplicate_count')!=pl.col('duplicate_count_source'))).height:raise ValueError('native source/catalog identity mismatch')
    total=sum(map(int,clones['duplicate_count']))
    if total>2**63-1:raise ValueError('native support count total exceeds safe integer bound')
    if entries['source_observation_id'].n_unique()!=entries.height:raise ValueError('one search segment per native catalog row required')
    if entries.join(catalog.select('source_observation_id'),on='source_observation_id',how='anti').height:raise ValueError('search entry outside native catalog')
    fields=['source_observation_id','native_clone_row','parent_clone_id','locus','native_junction_aa','duplicate_count']
    if entries.height:
        expected_rows=catalog.join(entries.select('source_observation_id'),on='source_observation_id',how='semi').select(fields).sort('source_observation_id')
        if not entries.select(fields).sort('source_observation_id').equals(expected_rows):raise ValueError('search entry provenance mismatch')
    records=records.join(entries.select('source_observation_id','seq'),on='source_observation_id',how='left',validate='1:1').with_columns(pl.col('seq').is_not_null().alias('search_entry_present'))
    missing=records.filter(~pl.col('search_entry_present'))
    if missing.filter((pl.col('native_junction_aa').str.len_chars()>=6)|(pl.col('coverage_kind')!='junction_only')).height:raise ValueError('missing search entry is only valid for short junction-only source rows')
    records=records.with_columns(pl.col('seq').fill_null(pl.col('native_junction_aa')).alias('seq'))
    if records.filter(~pl.col('seq').str.contains('^[ACDEFGHIKLMNPQRSTVWY]+$')|~pl.col('native_junction_aa').str.contains('^C[ACDEFGHIKLMNPQRSTVWY]+[FW]$')|(pl.col('seq').str.count_matches(pl.col('native_junction_aa'),literal=True)!=1)).height:raise ValueError('native sequence/junction unavailable or inconsistent')
    records=records.with_columns(pl.col('seq').str.find(pl.col('native_junction_aa'),literal=True).cast(pl.Int64).alias('junction_start')).with_columns((pl.col('junction_start')+pl.col('native_junction_aa').str.len_chars()).alias('junction_end'))
    records=records.select('source_observation_id','native_clone_row','parent_clone_id','locus','seq','native_junction_aa','junction_start','junction_end','coverage_kind','sequence_origin','background_member','search_entry_present','v_call','j_call','c_call',pl.col('duplicate_count').alias('native_support'))
    inventory=clones.group_by('locus').agg(pl.len().alias('source_native_clone_records'),pl.col('duplicate_count').sum().alias('source_native_support')).join(records.group_by('locus').agg(pl.len().alias('retained_catalog_records'),pl.col('native_support').sum().alias('retained_catalog_native_support')),on='locus',how='left').with_columns(pl.col('retained_catalog_records','retained_catalog_native_support').fill_null(0)).with_columns((pl.col('source_native_clone_records')-pl.col('retained_catalog_records')).alias('excluded_native_clone_records'))
    return receipt,records,inventory


def acquisition_calls(paths,geometry,observed_receipt):
    """Read exactly three native acquisitions once; share chymotrypsin through pool aliases."""
    import polars as pl
    from vdjaqp.readout import read_idxml
    blocks={};inventory=[];pins={};observed_hash=digest(observed_receipt);runs=set();physical_h=None
    for path in paths:
        path=Path(path);m=verify_manifest(path);block=m.get('block')
        if m.get('stage')!='acquisition' or m.get('geometry')!=geometry or block not in BLOCKS or block in blocks:raise ValueError('exactly three distinct same-geometry acquisition blocks required')
        if m['run_id'] in runs:raise ValueError('acquisition run IDs must be distinct')
        runs.add(m['run_id'])
        if not any(Path(p).name=='receipt.json' and h==observed_hash for p,h in m['input_sha256'].items()):raise ValueError('acquisition/observed source receipt differs')
        if physical_h is None:physical_h=m['physical_H']
        if m['physical_H']!=physical_h or physical_h<1 or (m['physical_L']<1 if block==BLOCKS[0] else m['physical_L']!=0):raise ValueError('physical preparation inventory differs')
        p=path.parent/'native/psms.idXML.gz'
        if m['artifact_sha256'].get('native/psms.idXML.gz')!=digest(p):raise ValueError('native PSM artifact missing or changed')
        with gzip.open(p,'rb') as f:calls=read_idxml(f,run_id=m['run_id'])
        blocks[block]=calls.select('peptide_id','spectrum_id','seq').with_columns(pl.lit(block).alias('acquisition'))
        pins.update({str(path.resolve()):digest(path),str(p.resolve()):digest(p)})
        inventory.append(dict(acquisition=block,run_id=m['run_id'],physical_H=m['physical_H'],physical_L=m['physical_L'],psm_assignments=calls.height,called_spectra=calls['spectrum_id'].n_unique(),distinct_called_AA=calls['seq'].n_unique(),used_in_solution=block!=BLOCKS[1],used_in_heavy_band=block!=BLOCKS[0]))
    if set(blocks)!=set(BLOCKS):raise ValueError('exactly three acquisition blocks required')
    calls=pl.concat([blocks[b] for b in BLOCKS])
    pools=pl.DataFrame({'acquisition':[BLOCKS[0],BLOCKS[2],BLOCKS[1],BLOCKS[2]],'condition':['solution','solution','heavy_band','heavy_band']})
    pooled=calls.join(pools,on='acquisition',validate='m:m')
    return pooled,pl.DataFrame(inventory),pins


def summarize(coverage,*,class_view=False):
    """Record and class-union weights have explicit, different support interpretations."""
    import polars as pl
    return coverage.group_by('condition','locus','region').agg(
        pl.len().alias('records'),pl.col('coverage').count().alias('defined_records'),
        (pl.col('covered_residues')>0).sum().alias('captured_records'),
        ((pl.col('covered_residues')>0)&(pl.col('covered_residues')<pl.col('region_length'))).sum().alias('partial_only_records'),
        ((pl.col('region_length')>0)&(pl.col('covered_residues')==pl.col('region_length'))).sum().alias('complete_records'),
        pl.col('covered_residues').sum().alias('covered_residues'),pl.col('region_length').sum().alias('defined_residues'),
        pl.col('native_support').sum().alias('native_support_denominator'),
        pl.col('native_support').filter(pl.col('covered_residues')>0).sum().alias('captured_native_support'),
        pl.col('native_support').filter(pl.col('coverage').is_not_null()).sum().alias('defined_native_support'),
        (pl.col('native_support')*pl.col('coverage')).sum().alias('_weighted')).with_columns(
            (pl.col('captured_records')/pl.col('records')).alias('captured_record_fraction'),
            (pl.col('complete_records')/pl.col('records')).alias('complete_record_fraction'),
            pl.when(pl.col('native_support_denominator')>0).then(pl.col('captured_native_support')/pl.col('native_support_denominator')).otherwise(None).alias('native_support_capture_fraction'),
            pl.when(pl.col('defined_native_support')>0).then(pl.col('_weighted')/pl.col('defined_native_support')).otherwise(None).alias('native_support_weighted_mean_coverage'),
            pl.lit(UNIT).alias('weight_unit'),pl.lit('member support of an exact-AA class with union capture; not each member independently captured' if class_view else 'native corrected clone record').alias('weight_scope')).drop('_weighted').sort('condition','locus','region')


def gene_tables(records,coverage):
    import polars as pl
    associated=coverage.filter(pl.col('region')=='OBSERVED_SEGMENT').select('source_observation_id','condition',(pl.col('covered_residues')>0).alias('_supported'))
    g=records.join(associated,on='source_observation_id',validate='1:m').unpivot(on=['v_call','j_call','c_call'],index=['source_observation_id','locus','native_support','condition','_supported'],variable_name='segment',value_name='native_call').with_columns(pl.col('segment').str.slice(0,1).str.to_uppercase())
    call=pl.col('native_call').fill_null('').str.strip_chars()
    missing=call.str.to_lowercase().is_in(['','none','null','nan','unknown','unassigned'])
    constant_prefix=pl.when(pl.col('locus')=='IGH').then(call.str.slice(3,1).is_in(['A','G','M','D','E'])).otherwise(call.str.slice(3,1)=='C')
    prefix=pl.when(pl.col('segment')=='C').then(call.str.starts_with(pl.col('locus')) & constant_prefix).otherwise(call.str.starts_with(pl.col('locus')+pl.col('segment')))
    single=call.str.contains(r'^IG[HKL][A-Za-z0-9/-]+(?:\*[A-Za-z0-9]+)?$') & prefix & (call.str.count_matches(r'IG[HKL]')==1)
    g=g.with_columns(pl.when(missing).then(pl.lit('missing')).when(single).then(pl.lit('single_call')).otherwise(pl.lit('compound_or_unresolved')).alias('call_status'),pl.when(missing).then(pl.lit('<missing>')).when(single).then(call.str.split('*').list.first()).otherwise(call).alias('call_bin'))
    bins=g.group_by('condition','locus','segment','call_bin','call_status').agg(pl.len().alias('native_records'),pl.col('_supported').sum().alias('associated_records'),pl.col('native_support').sum().alias('native_support_denominator'),pl.col('native_support').filter(pl.col('_supported')).sum().alias('associated_native_support')).with_columns((pl.col('associated_records')/pl.col('native_records')).alias('associated_record_fraction'),pl.when(pl.col('native_support_denominator')>0).then(pl.col('associated_native_support')/pl.col('native_support_denominator')).otherwise(None).alias('native_support_capture_fraction'),pl.lit(UNIT).alias('weight_unit'))
    summary=bins.group_by('condition','locus','segment').agg((pl.col('call_status')!='missing').sum().alias('distinct_called_bins'),((pl.col('call_status')!='missing')&(pl.col('associated_records')>0)).sum().alias('associated_captured_bins'),(pl.col('call_status')=='single_call').sum().alias('distinct_single_call_bins'),((pl.col('call_status')=='single_call')&(pl.col('associated_records')>0)).sum().alias('associated_captured_single_call_bins'),(pl.col('call_status')=='compound_or_unresolved').sum().alias('compound_or_unresolved_bins'),pl.col('native_records').filter(pl.col('call_status')=='missing').sum().alias('missing_call_records')).with_columns(pl.when(pl.col('distinct_called_bins')>0).then(pl.col('associated_captured_bins')/pl.col('distinct_called_bins')).otherwise(None).alias('associated_called_bin_fraction'),pl.when(pl.col('distinct_single_call_bins')>0).then(pl.col('associated_captured_single_call_bins')/pl.col('distinct_single_call_bins')).otherwise(None).alias('associated_single_call_fraction'))
    return bins.with_columns(pl.when(pl.col('segment')=='C').then(pl.lit('native_C_call_not_gene_or_subclass_verified')).otherwise(pl.lit('native_VJ_gene_call')).alias('call_namespace')).sort('condition','locus','segment','call_bin'),summary.sort('condition','locus','segment')


def capture(records,pooled_calls,*,cpus,max_placements,max_expanded_slots):
    """One native search; all remaining selection and residue unions are keyed tables."""
    import polars as pl
    from vdjaqp.simulation_metrics import exact_peptide_placements
    queries=pooled_calls.select('seq').unique().sort('seq').filter(pl.col('seq').str.len_chars()>=6).with_row_index('_q').with_columns(pl.format('peptide_{}',pl.col('_q')).alias('peptide_id')).drop('_q')
    excluded=pooled_calls.filter(pl.col('seq').str.len_chars()<6).select('acquisition','seq').unique()
    condition_queries=pooled_calls.select('seq','condition','spectrum_id').join(queries,on='seq').group_by('condition','peptide_id').agg(pl.len().alias('psm_assignments'),pl.col('spectrum_id').n_unique().alias('called_spectra'))
    targets=records.select(pl.col('seq').alias('protein_aa')).unique().sort('protein_aa').with_row_index('_t').with_columns(pl.format('target_{}',pl.col('_t')).alias('chain_id')).drop('_t')
    aliases=records.join(targets,left_on='seq',right_on='protein_aa',validate='m:1')
    hits=exact_peptide_placements(targets.select('chain_id','protein_aa'),queries,threads=cpus)
    if hits.height>max_placements:raise ValueError('native hit count exceeds explicit bound')
    estimate=(hits.join(aliases.group_by('chain_id').agg(pl.len().cast(pl.Int64).alias('_aliases')),on='chain_id').join(condition_queries.group_by('peptide_id').agg(pl.len().cast(pl.Int64).alias('_conditions')),on='peptide_id').select((pl.col('_aliases')*pl.col('_conditions')).sum()).item() or 0)
    if estimate>max_placements:raise ValueError('alias/condition association expansion exceeds explicit bound')
    placed=hits.join(aliases.select('chain_id','source_observation_id'),on='chain_id').join(condition_queries.select('condition','peptide_id'),on='peptide_id').select('source_observation_id','condition','peptide_id','start','end').unique()
    rg=pl.concat([records.select(pl.col('source_observation_id').alias('chain_id'),'locus',pl.lit(name).alias('region'),start.alias('start'),end.alias('end'),pl.lit('available_native_export').alias('status')) for name,start,end in [('OBSERVED_SEGMENT',pl.lit(0,dtype=pl.Int64),pl.col('seq').str.len_chars().cast(pl.Int64)),('JUNCTION',pl.col('junction_start'),pl.col('junction_end')),('CDR3',pl.col('junction_start')+1,pl.col('junction_end')-1)]])
    placements=placed.select(pl.col('source_observation_id').alias('chain_id'),pl.col('condition').alias('evidence_level'),'start','end').unique()
    coverage,cost=bounded_coverage(rg,placements,list(CONDITIONS),max_slots=max_expanded_slots)
    coverage=coverage.rename({'chain_id':'source_observation_id','evidence_level':'condition'}).join(records.select('source_observation_id','native_support','background_member','search_entry_present','coverage_kind'),on='source_observation_id',validate='m:1').with_columns(pl.lit(LEVEL).alias('evidence_level'))
    classes=records.select('locus','native_junction_aa').unique().sort('locus','native_junction_aa').with_row_index('_j').with_columns(pl.format('junction_{}',pl.col('_j')).alias('junction_class_id')).drop('_j')
    membership=records.join(classes,on=['locus','native_junction_aa'],validate='m:1')
    totals=membership.group_by('junction_class_id').agg(pl.len().alias('member_records'),pl.col('native_support').sum())
    classes=classes.join(totals,on='junction_class_id',validate='1:1')
    jc=placed.join(membership.select('source_observation_id','junction_class_id','junction_start','junction_end'),on='source_observation_id',validate='m:1').with_columns((pl.max_horizontal('start','junction_start')-pl.col('junction_start')).cast(pl.Int64).alias('_start'),(pl.min_horizontal('end','junction_end')-pl.col('junction_start')).cast(pl.Int64).alias('_end')).filter(pl.col('_end')>pl.col('_start')).select(pl.col('junction_class_id').alias('chain_id'),pl.col('condition').alias('evidence_level'),pl.col('_start').alias('start'),pl.col('_end').alias('end')).unique()
    cr=pl.concat([classes.select(pl.col('junction_class_id').alias('chain_id'),'locus',pl.lit(name).alias('region'),pl.lit(start,dtype=pl.Int64).alias('start'),(pl.col('native_junction_aa').str.len_chars().cast(pl.Int64)-trim).alias('end'),pl.lit('defined_exact_AA_class').alias('status')) for name,start,trim in [('JUNCTION',0,0),('CDR3',1,1)]])
    class_coverage,class_cost=bounded_coverage(cr,jc,list(CONDITIONS),max_slots=max_expanded_slots)
    class_coverage=class_coverage.rename({'chain_id':'junction_class_id','evidence_level':'condition'}).join(classes,on=['junction_class_id','locus'],validate='m:1').with_columns((pl.col('covered_residues')>0).alias('any_partial_capture'),((pl.col('covered_residues')>0)&(pl.col('covered_residues')<pl.col('region_length'))).alias('partial_only_capture'),(pl.col('covered_residues')==pl.col('region_length')).alias('complete_capture'),pl.lit('exact anchored AA sequence within locus; I/L-equivalent classes retained separately').alias('class_definition'),pl.lit(LEVEL).alias('evidence_level'))
    ambiguity=placed.join(membership.select('source_observation_id','locus','junction_class_id'),on='source_observation_id').group_by('condition','peptide_id').agg(pl.col('source_observation_id').n_unique().alias('compatible_native_records'),pl.col('junction_class_id').n_unique().alias('compatible_exact_junction_classes'),pl.col('locus').n_unique().alias('compatible_loci'))
    ambiguity=condition_queries.join(ambiguity,on=['condition','peptide_id'],how='left').with_columns(pl.col('compatible_native_records','compatible_exact_junction_classes','compatible_loci').fill_null(0))
    genes,gene_summary=gene_tables(records,coverage)
    return dict(records=records,queries=queries,excluded_short_queries=excluded,condition_queries=condition_queries,distinct_targets=targets,target_aliases=aliases.select('chain_id','source_observation_id','native_clone_row','locus','background_member'),native_exact_placements=hits,record_placements=placed,record_region_coverage=coverage,record_capture=summarize(coverage),junction_class_coverage=class_coverage,junction_class_capture=summarize(class_coverage,class_view=True),genes=genes,gene_capture=gene_summary,mapping_ambiguity=ambiguity,record_expansion_profile=cost,junction_expansion_profile=class_cost),dict(native_search_batches=int(bool(targets.height and queries.height)),native_hits=hits.height,alias_condition_association_upper_bound=estimate,max_placements=max_placements,max_expanded_slots=max_expanded_slots)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['observed-receipt','native-clones','output']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--geometry',choices=['amplicon','rnaseq'],required=True);p.add_argument('--acquisition-manifest',type=Path,action='append',required=True)
    p.add_argument('--cpus',type=int,default=16);p.add_argument('--max-placements',type=int,default=5000000);p.add_argument('--max-expanded-slots',type=int,default=2000000)
    a=p.parse_args(argv)
    if a.output.exists():raise FileExistsError('fresh output required')
    if a.max_placements<1 or a.max_expanded_slots<1:p.error('positive resource bounds required')
    from vdjaqp import resolve_cpu_budget
    cpus=resolve_cpu_budget(a.cpus);started=time.monotonic()
    import polars as pl
    pins={str(path.resolve()):digest(path) for path in [a.observed_receipt,a.native_clones,*a.acquisition_manifest]}
    receipt,records,inventory=load_observed(a.observed_receipt,a.native_clones,a.geometry)
    pins.update({str((a.observed_receipt.parent/name).resolve()):h for name,h in receipt['artifact_sha256'].items()})
    calls,acquisitions,acq_pins=acquisition_calls(a.acquisition_manifest,a.geometry,a.observed_receipt);pins.update(acq_pins)
    frames,limits=capture(records,calls,cpus=cpus,max_placements=a.max_placements,max_expanded_slots=a.max_expanded_slots)
    pre=receipt['preprocess']
    if any(type(pre[k]) is not int or pre[k]<0 for k in ['input_pairs','merged_pairs','unmerged_pairs','retained_output_records']):raise ValueError('preprocessing totals must be exact nonnegative integers')
    if pre['input_pairs']!=pre['merged_pairs']+pre['unmerged_pairs'] or pre['retained_output_records']!=pre['merged_pairs']+2*pre['unmerged_pairs']:raise ValueError('preprocessing unit/accounting mismatch')
    totals=pl.DataFrame({'quantity':['original_input_pairs','preprocessed_fragment_records','retained_catalog_native_support'],'count':[pre['input_pairs'],pre['retained_output_records'],sum(map(int,records['native_support']))],'unit':['raw_read_pairs','processed_fragment_records',UNIT]})
    frames.update(acquisitions=acquisitions,input_totals=totals,native_clone_inventory=inventory)
    if any(digest(Path(path))!=h for path,h in pins.items()):raise ValueError('report inputs changed')
    a.output.mkdir(parents=True)
    for name,frame in frames.items():frame.write_parquet(a.output/(name+'.parquet'),compression='zstd')
    code=[Path(__file__),Path(sys.modules['vdjaqp.simulation_metrics'].__file__),Path(sys.modules['vdjaqp.readout'].__file__),Path(sys.modules['scripts.export_mock_observation_metrics'].__file__),Path(sys.modules['scripts.run_matched_reference_ms'].__file__),Path(sys.modules['vdjaqp.match'].__file__)]
    report=dict(status='complete',geometry=a.geometry,cpu_budget=cpus,source_input_sha256=pins,observed_artifact_sha256=receipt['artifact_sha256'],observed_code_sha256=receipt.get('code_sha256',{}),observed_producer_provenance={k:receipt[k] for k in ['export_source_commit','native_producer_source_commit','native_producer_receipt_sha256'] if k in receipt},source_code_sha256={str(path.resolve()):digest(path) for path in code},rows={n:f.height for n,f in frames.items()},limits=limits,weight_unit=UNIT,denominator='valid native IGH/IGK/IGL catalog records only; inserted references excluded; native support is not original pair abundance',evidence_level=LEVEL,scope='associated native called bins (C is not verified gene/subclass) and observed segment/junction/CDR3 exact compatibility; no light scores, identity, direct V/J/C residue coverage or N-only inference',shared_acquisition='shared_chymotrypsin read once, contributes once per existing condition; conditions are not independent replicates',native_source_verification='supplied clone TSV byte-verified against canonical source pin; other native raw-source pins are historical producer provenance and not reverified here',artifact_sha256={f.name:digest(f) for f in a.output.iterdir()},wall_seconds=time.monotonic()-started)
    (a.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
