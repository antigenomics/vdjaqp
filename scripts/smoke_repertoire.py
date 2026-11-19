"""Extract observed receptor candidates with fastp and one native ARDA sample run.

Native AIRR projection/selection is lazy and vectorized; no generating truth or
unobserved framework is used. Original native TSVs retain unselected variants.
"""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

KEYS=['locus','v_call','j_call','junction']
# Raw corrected-clone geometry; never confuse these with read-space AIRR spans.
NATIVE_CLONE_GEOMETRY_FIELDS=('d_call','d2_call','d_support','d2_support','np1','np2','np3',
    'v_sequence_end','d_sequence_start','d_sequence_end','d2_sequence_start','d2_sequence_end','j_sequence_start')
NATIVE_CLONE_GEOMETRY_PROVENANCE='native inferred corrected-junction nucleotide geometry; 1-based closed coordinates; -1 means not located; raw strings, null when absent'
REGIONS=[('fwr1',0),('cdr1',1),('fwr2',2),('cdr2',3),('fwr3',4),('cdr3',5),('fwr4',7)]


def digest(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def _translate(expr):
    from arda.refbuild.translate import CODON_TABLE
    import polars as pl
    return expr.str.extract_all('...').list.eval(pl.element().replace_strict(CODON_TABLE,default='X')).list.join('')


def _scan(path, columns, index):
    import polars as pl
    frame=pl.scan_csv(path,separator='\t',quote_char=None,infer_schema_length=0,row_index_name=index)
    present=frame.collect_schema().names()
    return frame.select(pl.col(index).cast(pl.Int64),*[pl.when(pl.col(c)=='""').then(None).otherwise(pl.col(c)).alias(c) if c in present else pl.lit(None,dtype=pl.String).alias(c) for c in columns])


def _observations(path, assembled, membership):
    """Project required AIRR columns and join explicit native read membership."""
    import polars as pl
    from arda.rnaseq.correct import _COMPLETE
    columns=KEYS+['junction_aa','sequence_id','sequence','rev_comp','c_call','mmseqs2_score','mmseqs2_tend','mmseqs2_t_vjend','mmseqs2_qend']+[f'{r}_{edge}' for r,_ in REGIONS for edge in ('start','end')]
    frame=_scan(path,columns,'native_observation_row').with_columns(_COMPLETE.fill_null(False).alias('native_complete_junction'))
    # Native corrected-clone row identity closes the join; raw keys may change.
    frame=frame.join(membership.lazy(),on='sequence_id',how='inner')
    frame=frame.filter(pl.col('sequence').is_not_null() & pl.col('rev_comp').is_in(['T','F'])).with_columns(
        pl.when(pl.col('rev_comp')=='T').then(pl.col('sequence').str.to_uppercase().str.reverse().str.replace_many({'A':'T','T':'A','C':'G','G':'C'})).otherwise(pl.col('sequence').str.to_uppercase()).alias('_nt'),
        pl.lit(int(assembled)).alias('_assembled'),pl.lit('arda_native_assembled_read' if assembled else 'arda_observed_read').alias('sequence_origin'))
    frame=frame.filter(pl.col('_nt').str.count_matches(pl.col('junction'),literal=True)==1).with_columns(pl.col('_nt').str.find(pl.col('junction'),literal=True).cast(pl.Int64).alias('_jp'))
    frame=frame.with_columns((pl.col('_jp')%3).alias('_frame')).with_columns(_translate(pl.col('_nt').str.slice(pl.col('_frame'))).alias('_aa'),((pl.col('_jp')-pl.col('_frame'))//3).alias('_js'))
    # The unique observed nucleotide junction fixes its translated coordinate;
    # a synonymous copy elsewhere in the protein does not make it ambiguous.
    frame=frame.filter(pl.col('_aa').str.slice(pl.col('_js'),pl.col('junction_aa').str.len_chars())==pl.col('junction_aa'))
    intervals=[]
    def interval(lo,hi,label,order):
        return pl.when((lo>=0)&(hi>lo)&(hi<=pl.col('_aa').str.len_chars())).then(pl.struct(lo.alias('start'),hi.alias('end'),pl.lit(label).alias('label'),pl.lit(order).alias('order')))
    for order,(region,label) in enumerate(REGIONS):
        lo=((pl.col(region+'_start').cast(pl.Int64,strict=False)-1-pl.col('_frame')+2)//3).clip(lower_bound=0)
        hi=((pl.col(region+'_end').cast(pl.Int64,strict=False)-pl.col('_frame'))//3).clip(upper_bound=pl.col('_aa').str.len_chars())
        intervals.append(interval(lo,hi,label,order))
    const=(pl.col('c_call').is_not_null() & (pl.col('mmseqs2_tend').cast(pl.Float64,strict=False)>pl.col('mmseqs2_t_vjend').cast(pl.Float64,strict=False)) & (pl.col('mmseqs2_t_vjend').cast(pl.Float64,strict=False)>0))
    lo=(pl.col('fwr4_end').cast(pl.Int64,strict=False)-pl.col('_frame'))//3
    hi=((pl.col('mmseqs2_qend').cast(pl.Float64,strict=False).cast(pl.Int64)-pl.col('_frame'))//3).clip(upper_bound=pl.col('_aa').str.len_chars())
    intervals.append(pl.when(const).then(interval(lo,hi,8,7)))
    intervals.append(pl.when(pl.col('native_complete_junction')).then(interval(pl.col('_js')+1,pl.col('_js')+pl.col('junction_aa').str.len_chars()-1,5,8)))
    frame=frame.with_columns(pl.concat_list(intervals).list.eval(pl.element().filter(pl.element().is_not_null())).list.sort().alias('_intervals'))
    elem=pl.element()
    frame=frame.filter(pl.col('_intervals').list.len()>0).with_columns(
        pl.col('_intervals').list.first().struct.field('start').alias('_lo'),
        pl.col('_intervals').list.eval(elem.struct.field('end')).list.max().alias('_hi'),
        pl.col('_intervals').list.eval((elem.struct.field('start')>elem.struct.field('end').cum_max().shift(1)).fill_null(False)).list.any().alias('_gap'))
    frame=frame.filter(~pl.col('_gap') & (pl.col('_lo')<=pl.col('_js')) & (pl.col('_hi')>=pl.col('_js')+pl.col('junction_aa').str.len_chars())).with_columns(
        pl.col('_aa').str.slice(pl.col('_lo'),pl.col('_hi')-pl.col('_lo')).alias('seq'),
        pl.col('_nt').str.slice(pl.col('_frame')+3*pl.col('_lo'),3*(pl.col('_hi')-pl.col('_lo'))).alias('nt_seq'),
        (pl.col('_hi')-pl.col('_lo')).alias('_length'),pl.col('mmseqs2_score').cast(pl.Float64,strict=False).fill_null(0).alias('_score'))
    return frame.filter(pl.col('seq').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$'))



def junction_protection(row):
    """Protect inferred V/J-excluded core; unresolved geometry protects the interior."""
    from arda.rnaseq.correct import _CANONICAL_AA
    junction=row['junction_aa'];nt=row['junction'];locus=row['locus']
    if not isinstance(junction,str) or not re.fullmatch(_CANONICAL_AA,junction) or not isinstance(nt,str) or len(nt)%3:
        raise ValueError('junction protection requires a complete in-frame native junction with conserved anchors')
    labels=[4]+[5]*(len(junction)-2)+[7]
    calls=all(isinstance(row.get(call),str) and re.fullmatch(re.escape(locus)+segment+r'[0-9][A-Za-z0-9/-]*(?:\*[0-9]+)?',row[call]) for call,segment in [('v_call','V'),('j_call','J')])
    raw=[row.get(k) for k in ('v_sequence_end','j_sequence_start')]
    bounds=[int(x) if isinstance(x,str) and re.fullmatch(r'-?[0-9]+',x) else None for x in raw]
    ve,js=bounds
    resolved=calls and ve is not None and js is not None and 3<=ve<=js-1<=len(nt)-3
    left,right=(ve,js-1) if resolved else (None,None)
    for i in range(1,len(junction)-1):
        if not resolved or (left<right and 3*i<right and 3*i+3>left):labels[i]=6
    return labels,dict(junction_core_status='inferred_vj_excluded_core' if resolved else 'conservative_whole_junction_unresolved',
        junction_core_start_nt=left,junction_core_end_nt=right,
        junction_core_protection_provenance='region6 is a V/J reversion and background protection mask, including D; inferred corrected-junction 0-based half-open core, any-overlap codons; unresolved protects anchor-excluded interior; not N-only support')


def export_candidates(prefix, output, *, geometry='amplicon', raw_airr=None,
                      representative_policy='longest', all_loci=False):
    """Select observed native clone-keyed frames without constructing sequence.

    ``native_root`` optionally prefers an observed complete junction matching
    the corrected parent's NT junction, before the unchanged length/score ties.
    No match retains the original selection with explicit unresolved metadata.
    ``longest`` preserves the original default output schema and selection.
    ``all_loci`` additionally retains typed IGH/IGK/IGL candidate records in
    ``all_locus_candidates.parquet``; existing heavy-only outputs are unchanged.
    """
    import polars as pl
    if representative_policy not in ('longest','native_root'):
        raise ValueError('representative_policy must be longest or native_root')
    prefix,output=Path(prefix),Path(output)
    legacy=output.suffix=='.parquet';directory=output.with_suffix('') if legacy else output
    if output.exists() or directory.exists():raise FileExistsError('observed output already exists')
    sources=[Path(str(prefix)+'.clones.tsv'),Path(raw_airr) if raw_airr is not None else Path(str(prefix)+'.airr.tsv')]
    assembled=Path(str(prefix)+'.assembled.airr.tsv')
    if assembled.exists():sources.append(assembled)
    read_map=Path(str(prefix)+'.read_ids.tsv')
    if not read_map.is_file():raise FileNotFoundError('native .read_ids.tsv membership is required; activate the pinned ARDA overlay')
    sources.append(read_map)
    sources.extend(path for suffix in ['.arda.json','.stats.json','.stats.tsv'] if (path:=Path(str(prefix)+suffix)).exists())
    pins={str(p.resolve()):digest(p) for p in sources}
    native_clones=_scan(sources[0],KEYS+['junction_aa','c_call','duplicate_count',*NATIVE_CLONE_GEOMETRY_FIELDS],'native_clone_row').collect(engine='streaming')
    clones=native_clones.filter(pl.col('locus').is_in(['IGH','IGK','IGL']))
    clones=clones.with_columns(pl.col('duplicate_count').cast(pl.Int64,strict=True),_translate(pl.col('junction')).alias('_junction_translation'))
    if clones.filter(pl.col('duplicate_count').is_null() | (pl.col('duplicate_count')<0)).height:raise ValueError('native clone read counts must be nonnegative integers')
    valid=clones.filter(pl.col('junction').str.contains(r'^[ACGT]+$') & pl.col('junction_aa').str.contains(r'^C[ACDEFGHIKLMNPQRSTVWY]+[FW]$') & (pl.col('_junction_translation')==pl.col('junction_aa')) & (pl.col('junction').str.len_chars()%3==0))
    membership=pl.read_csv(read_map,separator='\t',infer_schema_length=0).select(
        'sequence_id',pl.col('clone_row').cast(pl.Int64,strict=True).alias('native_clone_row')).unique()
    if membership.null_count().sum_horizontal().item() or membership.join(native_clones.select('native_clone_row'),on='native_clone_row',how='anti').height:
        raise ValueError('native read map contains null/unknown clone membership')
    membership=membership.join(valid.select('native_clone_row'),on='native_clone_row',how='semi')
    member_links={};source_ids=[]
    for name,path in [('raw',sources[1]),*([('assembled',assembled)] if assembled.exists() else [])]:
        ids=_scan(path,['sequence_id'],'_source_row').select('sequence_id').collect(engine='streaming')
        source_ids.append(ids)
        links=ids.join(membership,on='sequence_id',how='inner')
        n_links,n_ids=links.select(pl.len(),pl.col('sequence_id').n_unique()).row(0)
        member_links['n_'+name+'_members_linked']=n_links
        member_links['n_'+name+'_read_ids_linked']=n_ids
    if membership.join(pl.concat(source_ids).unique(),on='sequence_id',how='anti').height:
        raise ValueError('native read map IDs are absent from raw/assembled source namespaces')
    frames=[_observations(sources[1],False,membership)]
    if assembled.exists():frames.append(_observations(assembled,True,membership))
    observed=pl.concat(frames,how='vertical_relaxed')
    groups=['native_clone_row']
    order=['native_complete_junction','_length','_score','sequence_id','native_observation_row','_assembled']
    reverse=[False,False,False,True,True,False]
    if representative_policy=='native_root':
        observed=observed.join(valid.lazy().select('native_clone_row',
            pl.col('junction').alias('_native_root_junction')),on='native_clone_row',how='left',validate='m:1')
        observed=observed.with_columns((pl.col('native_complete_junction') &
            (pl.col('junction')==pl.col('_native_root_junction'))).fill_null(False).alias('_native_root_match'))
        order.insert(1,'_native_root_match');reverse.insert(1,False)
    fields=[c for c in observed.collect_schema().names() if c not in groups]
    # Native top-k state retains one ranked representative per key, not Python read rows.
    summaries=[pl.struct(fields).top_k_by(order,k=1,reverse=reverse).first().alias('_best'),
        pl.len().alias('native_matching_observations'),pl.col('_nt').n_unique().alias('native_sequence_variants'),
        pl.col('native_complete_junction').sum().alias('native_complete_matching_observations')]
    if representative_policy=='native_root':
        summaries.extend([pl.col('_native_root_match').sum().alias('native_root_matching_observations'),
            pl.col('_length').filter(pl.col('native_complete_junction')).max().alias('_default_complete_length')])
    best=observed.group_by(groups).agg(summaries).unnest('_best').collect(engine='streaming')
    joined=valid.join(best,on=groups,how='left',suffix='_observation').sort('native_clone_row')
    catalog=[];candidates=[];entries=[];all_candidates=[]
    for row in joined.iter_rows(named=True):
        parent=f'{geometry}_clone_{row["native_clone_row"]}'
        observed=row.get('seq') is not None and bool(row.get('native_complete_junction'))
        seq=row['seq'] if observed else row['junction_aa'];nt=row['nt_seq'] if observed else row['junction']
        observed_junction=row.get('junction_observation',row['junction']) if observed else None
        observed_aa=row.get('junction_aa_observation',row['junction_aa']) if observed else None
        sequence_junction=observed_aa if observed else row['junction_aa']
        if observed:
            lo=row['_lo'];labels=[-1]*len(seq)
            for interval in sorted(row['_intervals'],key=lambda x:x['order']):
                for pos in range(max(lo,interval['start']),min(row['_hi'],interval['end'])):labels[pos-lo]=interval['label']
            if -1 in labels:raise ValueError('native projected frame has an unexplained gap')
            start=row['_js']-lo;end=start+len(sequence_junction);origin=row['sequence_origin'];kind='observed_variable_fragment'
        else:labels=[4]+[5]*(len(seq)-2)+[7];start=0;end=len(seq);origin='native_corrected_junction' if legacy else 'arda_corrected_junction';kind='junction_only'
        protection_row=dict(row,junction_aa=sequence_junction,junction=observed_junction if observed else row['junction'])
        if observed and observed_junction!=row['junction']:
            # Corrected-core coordinates must not silently annotate a different raw core.
            protection_row.update(v_sequence_end=None,j_sequence_start=None)
        junction_labels,protection=junction_protection(protection_row)
        if observed and observed_junction!=row['junction']:
            protection['junction_core_status']='observed_junction_differs_from_corrected_root'
        labels[start:end]=junction_labels
        observation=parent+('_assembled_' if row.get('_assembled') else '_read_')+str(row['native_observation_row']) if observed else parent+'_junction'
        common=dict(source_observation_id=observation,parent_clone_id=parent,locus=row['locus'],native_clone_row=row['native_clone_row'],duplicate_count=row['duplicate_count'],native_junction_nt=row['junction'],native_junction_aa=row['junction_aa'],sequence_origin=origin,coverage_kind=kind,background_member=True,observed_read_id=row.get('sequence_id'),native_matching_observations=row.get('native_matching_observations') or 0,native_complete_matching_observations=row.get('native_complete_matching_observations') or 0,native_sequence_variants=row.get('native_sequence_variants') or 0,observed_junction_nt=observed_junction,observed_junction_aa=observed_aa,native_membership_source='native_correct_airr_read_map',observed_locus=row.get('locus_observation',row['locus']) if observed else None,observed_v_call=row.get('v_call_observation',row['v_call']) if observed else None,observed_j_call=row.get('j_call_observation',row['j_call']) if observed else None,observed_c_call=row.get('c_call_observation',row['c_call']) if observed else None,coding_nt_start=row['_frame']+3*row['_lo'] if observed else None,coding_nt_end=row['_frame']+3*row['_hi'] if observed else None)
        if representative_policy=='native_root':
            root_match=bool(row.get('_native_root_match')) if observed else None
            common.update(representative_policy=representative_policy,
                representative_native_root_match=root_match,
                native_root_matching_observations=row.get('native_root_matching_observations') or 0,
                representative_selection_status=('observed_native_root' if root_match else
                    'observed_fallback_no_root_match' if observed else 'junction_only_no_complete_observation'))
        common.update(protection)
        common.update({'native_clone_'+field:row.get(field) for field in NATIVE_CLONE_GEOMETRY_FIELDS})
        common['native_clone_geometry_provenance']=NATIVE_CLONE_GEOMETRY_PROVENANCE
        candidate=dict(clone_id=parent,seq=seq,nt_seq=nt,cdr3_seq=sequence_junction,cdr3_start=start,cdr3_end=end,region_ids=labels,v_name=row['v_call'],j_name=row['j_call'],c_name=row['c_call'] or 'unknown',**common)
        catalog.append(dict(status='observed_frame' if observed else 'native_corrected_junction_only',**common))
        if row['locus']=='IGH':candidates.append(candidate)
        if all_loci:all_candidates.append(candidate)
        search_eligible=len(seq)>=6
        if representative_policy=='native_root':
            # Preserve the original parent-level search population even when
            # the preferred observed frame is too short to match a peptide.
            search_eligible=((row.get('_default_complete_length') or 0)>=6 if
                row.get('native_complete_matching_observations') else len(row['junction_aa'])>=6)
        if search_eligible:entries.append(dict(entry_id=parent+'_segment0',seq=seq,nt_seq=nt,**common))
    if not candidates:raise ValueError('No usable native heavy sequence; no framework reconstruction attempted')
    directory.mkdir(parents=True)
    geometry_schema={**{'native_clone_'+field:pl.String for field in NATIVE_CLONE_GEOMETRY_FIELDS},'junction_core_start_nt':pl.Int64,'junction_core_end_nt':pl.Int64}
    if representative_policy=='native_root':
        geometry_schema.update(representative_native_root_match=pl.Boolean,
            native_root_matching_observations=pl.Int64,representative_policy=pl.String,
            representative_selection_status=pl.String)
    tables={'catalog':pl.DataFrame(catalog,schema_overrides=geometry_schema,infer_schema_length=None),'candidates':pl.DataFrame(candidates,schema_overrides=geometry_schema,infer_schema_length=None),'search_entries':pl.DataFrame(entries,infer_schema_length=None,schema_overrides={**geometry_schema,'entry_id':pl.String,'seq':pl.String,'source_observation_id':pl.String,'parent_clone_id':pl.String,'locus':pl.String,'background_member':pl.Boolean})}
    if all_loci:
        tables['all_locus_candidates']=pl.DataFrame(all_candidates,schema_overrides=geometry_schema,infer_schema_length=None)
    candidate_path=output if legacy else directory/'candidates.parquet'
    for name,frame in tables.items():frame.write_parquet(candidate_path if name=='candidates' else directory/(name+'.parquet'),compression='zstd')
    if any(digest(Path(p))!=h for p,h in pins.items()):raise ValueError('native inputs changed during export')
    invalid=clones.filter(pl.col('locus')=='IGH').height-valid.filter(pl.col('locus')=='IGH').height
    report=dict(native_member_links=member_links,native_membership='read_map',status='complete',geometry=geometry,source_format='arda_fastq_native',counts={k:f.height for k,f in tables.items()},observed_locus_counts=dict(tables['catalog'].group_by('locus').len().iter_rows()),source_input_sha256=pins,code_sha256={path.name:digest(path) for path in [Path(__file__),Path(__file__).with_name('smoke_preprocess.py')]},artifact_sha256={p.name:digest(p) for p in directory.glob('*.parquet')},corrected_heavy_clones=clones.filter(pl.col('locus')=='IGH').height,retained_observed_candidates=len(candidates),missing_observed_sequence=invalid,invalid_native_junctions=invalid,junction_only_candidates=sum(r['coverage_kind']=='junction_only' for r in candidates),observed_variable_fragment_candidates=sum(r['coverage_kind']=='observed_variable_fragment' for r in candidates),candidate_aa_lengths=[len(r['seq']) for r in candidates],native_clones_without_observed_frame=sum(r['coverage_kind']=='junction_only' for r in catalog),heavy_clones_without_observed_frame=sum(r['coverage_kind']=='junction_only' for r in candidates),nucleotide_provenance='Native per-read sequences and annotations retain their producing-version provenance; raw FASTQ is retained separately. Coding offsets are zero-based half-open on the oriented native sequence. Corrected junction-only rows have no observed offsets.',native_stage3_representation='Per-member rescued junctions; current native assembled AIRR does not export full coding sequences or FR coordinates',selection='longest complete observed coding frame, native score, stable source ID/row; assembled status only breaks an otherwise exact tie; unselected variants retained in pinned native TSVs',background_membership='observed_source_no_inserted_references')
    if representative_policy=='native_root':
        report.update(representative_policy=representative_policy,
            representative_selection_counts=tables['catalog'].group_by('representative_selection_status').len().sort('representative_selection_status').to_dicts(),
            selection='complete observed frame matching corrected native NT junction first, then original longest/score/stable ties; no root match retains original selection and explicit unresolved metadata; no sequence replacement')
    if legacy:report['artifact_sha256'][candidate_path.name]=digest(candidate_path)
    (directory/'receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('scenario',choices=['amplicon','rnaseq']);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path);p.add_argument('--threads',type=int,default=None);p.add_argument('--stage-timeout',type=float,default=300.)
    p.add_argument('--representative-policy',choices=['longest','native_root'],default='longest',
        help='optional observed representative preference; native_root uses existing corrected-junction membership only')
    p.add_argument('--all-loci',action='store_true',help='also export typed heavy, kappa and lambda candidates; preserves heavy-only outputs')
    a=p.parse_args(argv)
    if not math.isfinite(a.stage_timeout) or a.stage_timeout<=0:p.error('stage-timeout must be finite and positive')
    import vdjaqp
    threads=vdjaqp.resolve_cpu_budget(a.threads);out=a.root/f'arda-{a.scenario}';output=a.output or a.root/f'{a.scenario}_observed'
    if out.exists() or output.exists():raise FileExistsError('stage output already exists')
    inputs=[a.root/f'{a.scenario}_{mate}.fastq.gz' for mate in ['R1','R2']]
    pins={str(f.resolve()):digest(f) for f in inputs};start=time.monotonic();out.mkdir(parents=True)
    from scripts.smoke_preprocess import preprocess_pairs
    merged=out/'observed_merged_and_unmerged.fastq.gz'
    preprocess=preprocess_pairs(*inputs,merged,threads=threads,fastp=Path(__file__).resolve().parents[1]/'.tool-envs/preprocess/bin/fastp',stage_timeout=a.stage_timeout)
    command=[str(Path(sys.executable).with_name('arda')),a.scenario,'--r1',str(merged),'--out-dir',str(out),'--out-prefix','smoke','--threads',str(threads)]
    subprocess.run(command,check=True,timeout=a.stage_timeout)
    if not Path(str(out/'smoke')+'.read_ids.tsv').is_file():raise FileNotFoundError('native membership map missing; activate pinned ARDA overlay')
    if any(digest(Path(f))!=h for f,h in pins.items()):raise ValueError('original paired FASTQ changed')
    report=export_candidates(out/'smoke',output,geometry=a.scenario,representative_policy=a.representative_policy,all_loci=a.all_loci)
    from importlib.metadata import version
    native_binaries=[Path(command[0]),Path(preprocess['command'][0])] if 'command' in preprocess else [Path(command[0])]
    native_binary_sha256={str(path.resolve()):digest(path) for path in native_binaries if path.is_file()}
    receipt=output/'receipt.json';report.update(native_package_version={'arda-mapper':version('arda-mapper')},native_binary_sha256=native_binary_sha256,merged_fastq_sha256={str(merged.resolve()):digest(merged)},merged_fastq_bytes=merged.stat().st_size,cpu_budget=threads,command=command,preprocess=preprocess,input_fastq_sha256=pins,wall_seconds=time.monotonic()-start)
    receipt.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
