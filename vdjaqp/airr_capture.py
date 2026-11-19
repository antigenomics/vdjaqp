"""Agreeing-residue AIRR capture exports; no inferred experimental read weights."""
import numpy as np
import polars as pl
from .alphabet import encode_frame
from .layout import expand_alignments
from .data import to_gene_level,REGION_NAMES

VIEWS=('all_compatible_alignments','positive_llr_alignments','identified_candidates_at20')


def capture_tables(evidence,scores,reference_ids,*,counts=None,count_unit=None,count_provenance=None,chunk_size=200000):
    c=to_gene_level(evidence.clonotypes);ids=set(c['clone_id']);refs=set(reference_ids)
    if not refs<=ids or c['clone_id'].n_unique()!=c.height or chunk_size<1:raise ValueError('invalid candidate/reference keys or chunk size')
    if set(scores['clone_id'])!=ids or scores['clone_id'].n_unique()!=c.height:raise ValueError('score population differs')
    weights={}
    placeholder=count_unit=='placeholder'
    if placeholder and counts is None:raise ValueError('placeholder mode requires a complete unit-count table')
    if counts is not None:
        if count_unit not in ('reads','UMIs','placeholder') or not isinstance(count_provenance,str) or not count_provenance.strip():raise ValueError('counts require explicit reads/UMIs/placeholder unit and provenance')
        if set(counts.columns)!={'clone_id','count'} or counts['clone_id'].n_unique()!=counts.height or set(counts['clone_id'])!=(ids if placeholder else ids-refs):raise ValueError('counts must cover the declared population: all IDs for placeholders, background IDs for measured counts')
        if counts['count'].dtype not in (pl.Int8,pl.Int16,pl.Int32,pl.Int64,pl.UInt8,pl.UInt16,pl.UInt32,pl.UInt64) or counts['count'].null_count():
            raise ValueError('counts require exact nonmissing integer storage, not floating point')
        v=counts['count'].to_numpy()
        if placeholder and not (v==1).all():raise ValueError('placeholder counts must all equal one')
        total=sum(map(int,v))
        if (v<0).any() or not 0<total<=np.iinfo(np.int64).max:raise ValueError('counts negative, empty or overflow safe integer total')
        weights=dict(zip(counts['clone_id'],map(int,v)))
    code,offset,length=encode_frame(c['seq']);pc,po,plength=encode_frame(evidence.peptides['seq'])
    region=np.concatenate(c['region_ids'].to_list()).astype(int)
    if not np.array_equal(c['region_ids'].list.len().to_numpy(),length) or len(region)!=len(code) or not np.isin(region,list(REGION_NAMES)).all():raise ValueError('invalid supplied region labels')
    start=c['cdr3_start'].to_numpy();end=c['cdr3_end'].to_numpy()
    if c.filter((pl.col('cdr3_start')<0)|(pl.col('cdr3_end')<=pl.col('cdr3_start'))|
                (pl.col('cdr3_end')>pl.col('seq').str.len_chars())|
                (pl.col('seq').str.slice(pl.col('cdr3_start'),pl.col('cdr3_end')-pl.col('cdr3_start'))!=pl.col('cdr3_seq'))).height:
        raise ValueError('junction does not match supplied sequence slice')
    core_fields=['junction_core_status','junction_core_start_nt','junction_core_end_nt','junction_core_protection_provenance']
    native=(c['junction_core_status'].is_not_null().to_numpy() if 'junction_core_status' in c.columns
            else np.zeros(c.height,dtype=bool))
    has_native=bool(native.any())
    if has_native and not set(core_fields)<=set(c.columns):
        raise ValueError('native junction protection requires status, bounds and provenance')
    clone_row=np.repeat(np.arange(c.height),length);within=np.arange(len(code))-offset[clone_row]
    junction=(within>=start[clone_row])&(within<end[clone_row])
    keyed=evidence.alignments.join(c.select('clone_id').with_row_index('_c'),on='clone_id',validate='m:1').join(evidence.peptides.select('peptide_id').with_row_index('_p'),on='peptide_id',validate='m:1')
    selected=set(scores.filter(pl.col('confirmed_nats')>=20)['clone_id'])
    coverage=np.zeros((3,len(code)),dtype=bool)
    for begin in range(0,keyed.height,chunk_size):
        a=keyed.slice(begin,chunk_size);cr=a['_c'].to_numpy();pr=a['_p'].to_numpy();st=a['start'].to_numpy()
        if (st<0).any() or (st+plength[pr]>length[cr]).any():raise ValueError('alignment outside supplied sequence')
        ai,ps,rs=expand_alignments(po,plength,pr,offset,cr,st)
        agree=(pc[ps]>=0)&(pc[ps]==code[rs])
        for view,keep in enumerate([np.ones(a.height,bool),a['llr'].to_numpy()>0,a['clone_id'].is_in(selected).to_numpy()]):coverage[view,rs[agree & keep[ai]]]=True
    base=c.select('clone_id','cdr3_seq','v_name','j_name','c_name',*(core_fields if has_native else [])).with_columns(pl.Series('is_reference',[x in refs for x in c['clone_id']]),pl.Series('airr_count',[weights.get(x) for x in c['clone_id']],dtype=pl.UInt64))
    junction_classes=c.select(pl.col('cdr3_seq').alias('junction_aa')).unique(maintain_order=True).with_row_index('_j')
    jc=c.select(pl.col('cdr3_seq').alias('junction_aa')).join(junction_classes,on='junction_aa',how='left',validate='m:1')['_j'].to_numpy()
    jl=junction_classes['junction_aa'].str.len_bytes().to_numpy().astype(np.int64);jo=np.r_[0,np.cumsum(jl)]
    junction_slots=jo[jc[clone_row[junction]]]+within[junction]-start[clone_row[junction]]
    junction_owner=np.repeat(np.arange(len(jl)),jl)
    records=[];gene_rows=[];junction_rows=[];ambiguity=[]
    masks={'protein':np.ones(len(code),bool),**{name:region==k for k,name in REGION_NAMES.items()},'anchor_inclusive_junction':junction,'supplied_V_framework_proxy_0_to_4':region<=4,'supplied_J_FR4_proxy_7':region==7,'supplied_C_stub_proxy_8':region==8}
    if has_native:
        masks['CDR3']=np.where(native[clone_row],
            (within>start[clone_row])&(within<end[clone_row]-1),region==5)
    for vi,view in enumerate(VIEWS):
        selected_hits=keyed if vi==0 else keyed.filter(pl.col('llr')>0) if vi==1 else keyed.filter(pl.col('clone_id').is_in(selected))
        ambiguity.append(selected_hits.join(c.select('clone_id','cdr3_seq'),on='clone_id',validate='m:1').group_by('peptide_id').agg(pl.col('clone_id').n_unique().alias('compatible_candidate_records'),pl.col('cdr3_seq').n_unique().alias('compatible_anchored_junction_strings')).with_columns(pl.lit(view).alias('view')))
        per={}
        for name,mask in masks.items():
            den=np.bincount(clone_row[mask],minlength=c.height);num=np.bincount(clone_row[mask & coverage[vi]],minlength=c.height);per[name]=(num,den)
            name_expr=(pl.when(pl.col('junction_core_status').is_not_null()).then(pl.lit('junction_protection')).otherwise(pl.lit(name))
                       if has_native and name=='NDN' else pl.lit(name))
            records.append(base.with_columns(pl.lit(view).alias('view'),name_expr.alias('region'),pl.Series('covered_aa',num),pl.Series('defined_aa',den)).with_columns(pl.when(pl.col('defined_aa')>0).then(pl.col('covered_aa')/pl.col('defined_aa')).otherwise(None).alias('coverage')))
        # Junction-position union preserves shared partial support without unique identity claims.
        for population,keep in [('all',np.ones(c.height,bool)),('background',~base['is_reference'].to_numpy()),('reference',base['is_reference'].to_numpy())]:
            union=np.zeros(int(jo[-1]),bool)
            eligible=keep[clone_row[junction]] & coverage[vi,junction]
            np.logical_or.at(union,junction_slots[eligible],True)
            covered=np.bincount(junction_owner[union],minlength=len(jl))
            membership=base.filter(pl.Series(keep)).group_by('cdr3_seq').agg(pl.len().alias('member_records'),pl.col('airr_count').sum())
            classes=membership.rename({'cdr3_seq':'junction_aa'}).join(junction_classes,on='junction_aa',validate='1:1')
            indices=classes['_j'].to_numpy()
            classes=classes.with_columns(pl.Series('covered_aa',covered[indices]),pl.Series('defined_aa',jl[indices]),
                pl.lit(view).alias('view'),pl.lit(population).alias('population'),pl.lit('exact anchored AA string; not unique biological clonotype').alias('junction_definition'))
            if counts is None or (population!='background' and not placeholder):classes=classes.with_columns(pl.lit(None,dtype=pl.UInt64).alias('airr_count'))
            junction_rows.append(classes.drop('_j').with_columns((pl.col('covered_aa')>0).alias('any_partial_capture'),(pl.col('covered_aa')==pl.col('defined_aa')).alias('complete_capture')))
            for segment,column,mask in [('V','v_name','supplied_V_framework_proxy_0_to_4'),('J','j_name','supplied_J_FR4_proxy_7'),('C','c_name','supplied_C_stub_proxy_8')]:
                frame=base.with_columns(pl.Series('associated',per['protein'][0]>0),pl.Series('direct_supplied_region',per[mask][0]>0)).filter(pl.Series(keep)).with_columns(pl.col(column).fill_null('').alias('gene_call'))
                for row in frame.group_by('gene_call').agg(pl.len().alias('records'),pl.col('associated').sum().alias('associated_records'),pl.col('direct_supplied_region').sum().alias('direct_region_records'),pl.col('airr_count').sum().alias('_weights'),pl.col('airr_count').filter(pl.col('associated')).sum().alias('_captured_weights')).iter_rows(named=True):
                    gene_rows.append(dict(view=view,population=population,segment=segment,gene_call=row['gene_call'],missing_gene=row['gene_call'].strip().lower() in ('','nan','none','null'),records=row['records'],associated_records=row['associated_records'],direct_region_records=row['direct_region_records'],record_fraction=row['associated_records']/row['records'],airr_count_denominator=row['_weights'] if counts is not None and (population=='background' or placeholder) else None,airr_captured_count=row['_captured_weights'] if counts is not None and (population=='background' or placeholder) else None,airr_capture_fraction=row['_captured_weights']/row['_weights'] if counts is not None and (population=='background' or placeholder) and row['_weights'] else None))
    records=pl.concat(records)
    totals=[]
    for population in ('all','background','reference'):
        frame=records if population=='all' else records.filter(~pl.col('is_reference') if population=='background' else pl.col('is_reference'))
        groups=frame.partition_by(['view','region'],as_dict=True) if frame.height else {(v,r):frame for v in VIEWS for r in masks}
        for (view,region_name),group in groups.items():
            n=group.height;defined=group.filter(pl.col('defined_aa')>0);captured=int((group['covered_aa']>0).sum());weighted=counts is not None and (population=='background' or placeholder)
            mass=int(group['airr_count'].sum()) if weighted else None
            cmass=int(group.filter(pl.col('covered_aa')>0)['airr_count'].sum()) if weighted else None
            denom=int(defined['airr_count'].sum()) if weighted else None
            totals.append(dict(view=view,region=region_name,population=population,records=n,defined_records=defined.height,captured_records=captured,captured_record_fraction=captured/n if n else None,captured_defined_record_fraction=captured/defined.height if defined.height else None,
                covered_aa=int(group['covered_aa'].sum()),defined_aa=int(group['defined_aa'].sum()),airr_count=mass,captured_airr_count=cmass,airr_capture_fraction=cmass/mass if mass else None,
                weighted_mean_coverage=float((defined['airr_count']*defined['coverage']).sum())/denom if denom else None))
    genes=pl.DataFrame(gene_rows,schema_overrides={"airr_count_denominator":pl.UInt64,"airr_captured_count":pl.UInt64});summary=genes.group_by('view','population','segment').agg((~pl.col('missing_gene')).sum().alias('distinct_called_bins'),((pl.col('associated_records')>0)&~pl.col('missing_gene')).sum().alias('associated_captured_bins'),((pl.col('direct_region_records')>0)&~pl.col('missing_gene')).sum().alias('direct_region_captured_bins')).with_columns(pl.when(pl.col('distinct_called_bins')>0).then(pl.col('associated_captured_bins')/pl.col('distinct_called_bins')).otherwise(None).alias('associated_capture_fraction'))
    out=dict(mapping_ambiguity=pl.concat(ambiguity),records=records,record_capture=pl.DataFrame(totals,schema_overrides={"airr_count":pl.UInt64,"captured_airr_count":pl.UInt64}),junctions=pl.concat(junction_rows),genes=genes,gene_capture=summary)
    if placeholder:
        out={name:frame.with_columns(pl.lit('placeholder').alias('count_unit'),pl.lit(True).alias('counts_are_placeholder'),pl.lit(count_provenance).alias('count_provenance')) for name,frame in out.items()}
    return out
