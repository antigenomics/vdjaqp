"""Evaluation-only physical coverage and native contributor correctness.

No simulation truth enters inference. Coordinates are protein, zero-based half-open.
"""
from __future__ import annotations
import re
from collections.abc import Sequence
import polars as pl
from .match import PeptideMatcher, MIN_QUERY_LEN, ms_collapse

_AA=re.compile('[ACDEFGHIKLMNPQRSTVWY]+')


def _keys(frame, columns):
    if any(c not in frame.columns for c in columns):
        raise ValueError(f'missing required columns: {columns}')
    if frame.select(columns).null_count().sum_horizontal().item() or frame.select(columns).is_duplicated().any():
        raise ValueError(f'null or duplicate keys: {columns}')


def exact_peptide_placements(chains: pl.DataFrame, peptides: pl.DataFrame, threads=None) -> pl.DataFrame:
    """All exact I/L-compatible placements via one native batch; reject short queries."""
    _keys(chains,['chain_id']); _keys(peptides,['peptide_id'])
    for column,frame,minimum in [('protein_aa',chains,1),('seq',peptides,MIN_QUERY_LEN)]:
        if column not in frame.columns or any(not isinstance(s,str) or len(s)<minimum or not _AA.fullmatch(s) for s in frame[column]):
            raise ValueError(f'unsearchable {column}: canonical sequence and minimum length {minimum} required')
    schema={'chain_id':pl.String,'peptide_id':pl.String,'start':pl.Int32,'end':pl.Int32}
    if not chains.height or not peptides.height:return pl.DataFrame(schema=schema)
    return PeptideMatcher(chains['chain_id'],chains['protein_aa']).search(
        peptides['peptide_id'],peptides['seq'],max_subs=0,max_indels=0,threads=0 if threads is None else threads
    ).rename({'clone_id':'chain_id'}).select(list(schema)).sort(['chain_id','peptide_id','start'])


def region_coverage(regions: pl.DataFrame, placements: pl.DataFrame, evidence_levels: Sequence[str]) -> pl.DataFrame:
    """Residue union per supplied full region and evidence level, preserving metadata.

    Supply explicit rows with null start/end for undefined regions. Missing loci
    cannot be inferred from absent rows. Multiple named regions may overlap (junction).
    """
    if 'status' in regions.columns:
        if 'region_status' in regions.columns:raise ValueError('ambiguous region status columns')
        regions=regions.rename({'status':'region_status'})
    _keys(regions,['chain_id','region'])
    required=['chain_id','start','end','evidence_level']
    if any(c not in placements.columns for c in required) or any(c not in regions.columns for c in ['locus','start','end']):
        raise ValueError('missing coordinate/level columns')
    levels=list(evidence_levels)
    if not levels or len(set(levels))!=len(levels) or any(not isinstance(s,str) or not s for s in levels):raise ValueError('invalid evidence levels')
    for frame in [regions,placements]:
        for name in ['start','end']:
            if frame[name].dtype!=pl.Null and not frame[name].dtype.is_integer():raise ValueError('coordinates must be integers')
    bad=regions.filter((pl.col('start').is_null()!=pl.col('end').is_null())|(pl.col('start')<0)|(pl.col('end')<pl.col('start')))
    if bad.height:raise ValueError('invalid region intervals')
    if placements.select(required).null_count().sum_horizontal().item() or placements.filter((pl.col('start')<0)|(pl.col('end')<=pl.col('start'))|~pl.col('evidence_level').is_in(levels)).height:
        raise ValueError('invalid placement interval or evidence level')
    bounds=regions.group_by('chain_id').agg(pl.col('end').max().alias('_length'))
    check=placements.join(bounds,on='chain_id',how='left',validate='m:1')
    if check.filter(pl.col('_length').is_null()|(pl.col('end')>pl.col('_length'))).height:raise ValueError('unknown chain or out-of-range placement')
    # Native table expansion/union; no Python loop per alignment or residue.
    slots=(placements.select(required).with_columns(pl.int_ranges('start','end').alias('_position'))
           .explode('_position',empty_as_null=False).select('chain_id','evidence_level','_position').unique())
    base=regions.with_columns((pl.col('end')-pl.col('start')).alias('region_length')).join(pl.DataFrame({'evidence_level':levels}),how='cross')
    valid=(base.filter(pl.col('region_length')>0).with_columns(pl.int_ranges('start','end').alias('_position')).explode('_position',empty_as_null=False))
    counts=(valid.join(slots,on=['chain_id','evidence_level','_position'],how='inner')
            .group_by(['chain_id','region','evidence_level']).agg(pl.len().cast(pl.Int64).alias('covered_residues')))
    return (base.join(counts,on=['chain_id','region','evidence_level'],how='left',validate='1:1')
            .with_columns(pl.when(pl.col('region_length')>0).then(pl.col('covered_residues').fill_null(0)).otherwise(None).alias('covered_residues'))
            .with_columns((pl.col('covered_residues')/pl.col('region_length')).alias('coverage'),
                          pl.when(pl.col('region_length')>0).then(pl.lit('defined')).otherwise(pl.lit('undefined')).alias('status'))
            .sort(['chain_id','region','evidence_level']))


def psm_correctness(peptides: pl.DataFrame, native_truth: Sequence[dict]) -> tuple[pl.DataFrame,pl.DataFrame]:
    """PSM correctness against contributor sets plus separate assignment/event counts.

    Missing native links are unscorable (null correctness), never counted wrong.
    Matching protein count is compatibility, not unique antibody recovery.
    """
    _keys(peptides,['peptide_id'])
    if any(c not in peptides.columns for c in ['seq','native_spectrum_reference']):raise ValueError('missing native PSM fields')
    if any(not isinstance(s,str) or not _AA.fullmatch(s) for s in peptides['seq']):raise ValueError('noncanonical called peptide')
    if any(not isinstance(s,str) or not s for s in peptides['native_spectrum_reference']):raise ValueError('missing native spectrum reference')
    events=[];contributors=[];proteins=[];seen=set()
    for event in native_truth:
        key=event['native_spectrum_reference']
        if not isinstance(key,str) or not key or key in seen:raise ValueError('duplicate/missing native event identity')
        seen.add(key)
        items=event['contributors']
        events.append((key,len(items)))
        for item in items:
            seq=item['peptide_aa']
            if not isinstance(seq,str) or not _AA.fullmatch(seq):raise ValueError('noncanonical native contributor peptide')
            normalized=ms_collapse([seq])[0]
            contributors.append((key,normalized))
            proteins.extend((key,normalized,p) for p in set(item['protein_accessions']))
    event_frame=pl.DataFrame(events,schema={'native_spectrum_reference':pl.String,'native_contributor_count':pl.Int64},orient='row')
    truth=pl.DataFrame(contributors,schema={'native_spectrum_reference':pl.String,'_normalized':pl.String},orient='row').unique()
    protein_frame=pl.DataFrame(proteins,schema={'native_spectrum_reference':pl.String,'_normalized':pl.String,'protein':pl.String},orient='row').unique()
    matches=protein_frame.group_by(['native_spectrum_reference','_normalized']).agg(pl.len().cast(pl.Int64).alias('matching_protein_count'))
    annotated=(peptides.with_row_index('_psm_order').with_columns(pl.col('seq').str.replace_all('I','L',literal=True).alias('_normalized'))
        .join(event_frame,on='native_spectrum_reference',how='left',validate='m:1')
        .join(truth.with_columns(pl.lit(True).alias('_match')),on=['native_spectrum_reference','_normalized'],how='left',validate='m:1')
        .join(matches,on=['native_spectrum_reference','_normalized'],how='left',validate='m:1')
        .with_columns(pl.when(pl.col('native_contributor_count').is_not_null()).then(pl.col('_match').fill_null(False)).otherwise(None).alias('correct'),
                      pl.when(pl.col('native_contributor_count').is_not_null()).then(pl.col('matching_protein_count').fill_null(0)).otherwise(None).alias('matching_protein_count'))
        .with_columns(pl.when(pl.col('correct').is_null()).then(pl.lit('missing_native_truth')).otherwise(pl.lit('evaluated')).alias('status'))
        .with_columns((pl.col('matching_protein_count')>1).alias('protein_ambiguous'))
        .sort('_psm_order').drop('_normalized','_match','_psm_order'))
    evaluated=annotated.filter(pl.col('correct').is_not_null())
    counts={'psm_assignments':peptides.height,'native_ms2_events':len(events),'evaluated_psm_assignments':evaluated.height,
            'unscorable_psm_assignments':peptides.height-evaluated.height,'correct_psm_assignments':int(evaluated['correct'].sum() or 0),
            'called_spectra':peptides['native_spectrum_reference'].n_unique(),'evaluated_spectra':evaluated['native_spectrum_reference'].n_unique(),
            'spectra_with_any_correct_call':evaluated.filter(pl.col('correct'))['native_spectrum_reference'].n_unique(),
            'distinct_called_peptides':peptides['seq'].n_unique(),'distinct_il_equivalent_called_peptides':peptides['seq'].str.replace_all('I','L',literal=True).n_unique()}
    rows=[(name,count,None,None,'count') for name,count in counts.items()]
    den=evaluated.height;num=counts['correct_psm_assignments']
    rows.append(('correct_psm_rate',num,den,num/den if den else None,'fraction'))
    summary=pl.DataFrame(rows,schema={'quantity':pl.String,'count':pl.Int64,'denominator':pl.Int64,'fraction':pl.Float64,'unit':pl.String},orient='row')
    return annotated,summary


def identity_metrics(metadata: pl.DataFrame, candidates: pl.DataFrame,
                     ranks: pl.DataFrame, associations: pl.DataFrame):
    """Heavy-variable mass-class recovery at the existing 20 confirmed-nat criterion.

    Exact fragment associations are evaluation-only. A partial fragment never
    resolves a full variable identity. Physical rows and mass classes remain distinct.
    """
    import hashlib
    _keys(metadata,['chain_id']); _keys(candidates,['clone_id']); _keys(ranks,['clone_id'])
    if set(candidates.columns)&{'rank','n_tied','confirmed_nats'}:raise ValueError('candidate frame must not shadow authoritative rank columns')
    for name in ['rank','n_tied','confirmed_nats']:
        if name not in ranks.columns:raise ValueError('missing rank/evidence fields')
    if (ranks.select('clone_id').join(candidates.select('clone_id'),on='clone_id',how='anti').height or
        candidates.select('clone_id').join(ranks.select('clone_id'),on='clone_id',how='anti').height):
        raise ValueError('candidate and rank populations differ')
    if any(not ranks[c].dtype.is_integer() for c in ['rank','n_tied']) or ranks.filter(
        pl.col('rank').is_null()|pl.col('n_tied').is_null()|(pl.col('rank')<1)|(pl.col('rank')>ranks.height)|(pl.col('n_tied')<1)|(pl.col('n_tied')>ranks.height)|pl.col('confirmed_nats').is_null()|~pl.col('confirmed_nats').is_finite()).height:
        raise ValueError('invalid ranks, tie counts or confirmed evidence')
    if ranks['rank'].n_unique()!=ranks.height:raise ValueError('duplicate deterministic ranks')
    heavy=metadata.filter(pl.col('locus')=='IGH')
    if any(not isinstance(s,str) or not _AA.fullmatch(s) for s in heavy['variable_aa']):raise ValueError('invalid full variable sequence')
    if 'variable_length' not in heavy.columns:heavy=heavy.with_columns(pl.col('variable_aa').str.len_chars().alias('variable_length'))
    if heavy.filter(pl.col('variable_length').is_null()|(pl.col('variable_length')!=pl.col('variable_aa').str.len_chars())).height:raise ValueError('invalid variable lengths')
    labels=['mass_'+hashlib.sha256(s.encode()).hexdigest() for s in ms_collapse(heavy['variable_aa'])]
    heavy=heavy.with_columns(pl.Series('variable_mass_class_id',labels,dtype=pl.String))
    needed=['clone_id','chain_id','start','end']
    if any(c not in associations.columns for c in needed):raise ValueError('missing association columns')
    a=associations.select(needed).unique()
    if a.select(needed).null_count().sum_horizontal().item() or any(not a[c].dtype.is_integer() for c in ['start','end']):raise ValueError('invalid association coordinates')
    if a.join(candidates.select('clone_id'),on='clone_id',how='anti').height or a.join(heavy.select('chain_id'),on='chain_id',how='anti').height:raise ValueError('unknown candidate or non-heavy truth key')
    bounds=heavy.select('chain_id','variable_mass_class_id','variable_length',pl.col('protein_aa').str.len_chars().alias('_protein_length'))
    a=a.join(bounds,on='chain_id',how='left',validate='m:1')
    if a.filter((pl.col('start')<0)|(pl.col('end')<=pl.col('start'))|(pl.col('end')>pl.col('_protein_length'))).height:raise ValueError('association outside physical chain')
    pairs=(a.with_columns(((pl.col('start')==0)&(pl.col('end')>=pl.col('variable_length'))).alias('_full'))
        .group_by(['clone_id','variable_mass_class_id']).agg(pl.col('_full').any()))
    cc=pairs.group_by('clone_id').agg(pl.col('_full').sum().cast(pl.Int64).alias('full_class_count'),(~pl.col('_full')).sum().cast(pl.Int64).alias('partial_class_count'))
    called=(candidates.join(ranks.select('clone_id','rank','n_tied','confirmed_nats'),on='clone_id',how='left',validate='1:1')
        .join(cc,on='clone_id',how='left',validate='1:1').with_columns(pl.col('full_class_count','partial_class_count').fill_null(0),
            (pl.col('confirmed_nats')>=20).alias('shortlisted_at20'))
        .with_columns(pl.when(pl.col('full_class_count')==1).then(pl.lit('unique_full')).when(pl.col('full_class_count')>1).then(pl.lit('ambiguous_full')).when(pl.col('partial_class_count')>0).then(pl.lit('partial_only')).otherwise(pl.lit('no_exact_compatibility')).alias('identity_status')))
    full=pairs.filter(pl.col('_full')).join(called.select('clone_id','rank','n_tied','full_class_count','shortlisted_at20'),on='clone_id',validate='m:1')
    best=(full.sort(['rank','clone_id']).group_by('variable_mass_class_id',maintain_order=True)
          .agg(pl.col('clone_id').first().alias('best_full_compatible_candidate_id'),pl.col('rank').first().alias('best_full_compatible_rank'),pl.col('n_tied').first().alias('best_full_compatible_n_tied')))
    totals=pairs.group_by('variable_mass_class_id').agg(pl.col('_full').sum().cast(pl.Int64).alias('full_variable_compatible_candidate_count'),(~pl.col('_full')).sum().cast(pl.Int64).alias('partial_compatible_candidate_count'))
    recovered=full.group_by('variable_mass_class_id').agg(((pl.col('full_class_count')==1)&pl.col('shortlisted_at20')).any().alias('class_recovered_at20'))
    population=(heavy.join(totals,on='variable_mass_class_id',how='left',validate='m:1')
        .join(best,on='variable_mass_class_id',how='left',validate='m:1').join(recovered,on='variable_mass_class_id',how='left',validate='m:1')
        .with_columns(pl.col('full_variable_compatible_candidate_count','partial_compatible_candidate_count').fill_null(0),pl.col('class_recovered_at20').fill_null(False)))
    classes=population.unique('variable_mass_class_id')
    observable=classes.filter(pl.col('full_variable_compatible_candidate_count')>0).height
    n_recovered=classes.filter(pl.col('class_recovered_at20')).height
    shortlist=called.filter(pl.col('shortlisted_at20'))
    known=shortlist.filter(pl.col('full_class_count')==1).height
    counts={'true_heavy_rows':heavy.height,'true_variable_mass_classes':classes.height,'full_variable_observable_classes':observable,
            'shortlisted_candidates_at20':shortlist.height,'unique_full_true_candidates_in_shortlist':known,'unresolved_shortlist_candidates':shortlist.height-known}
    rows=[(k,v,None,None,'count') for k,v in counts.items()]
    rows.extend([('recovered_classes_at20',n_recovered,classes.height,n_recovered/classes.height if classes.height else None,'fraction'),
                 ('conditional_class_recovery_at20',n_recovered,observable,n_recovered/observable if observable else None,'fraction')])
    summary=pl.DataFrame(rows,schema={'quantity':pl.String,'count':pl.Int64,'denominator':pl.Int64,'fraction':pl.Float64,'unit':pl.String},orient='row')
    return population,called,summary
