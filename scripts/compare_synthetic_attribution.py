"""Fresh synthetic scoring comparison; final observations are inputs, truth is evaluation-only."""
from pathlib import Path
import csv,gzip,hashlib,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from vdjaqp import resolve_cpu_budget

CANDIDATE_COLUMNS=('clone_id','seq','nt_seq','cdr3_seq','cdr3_start','cdr3_end','region_ids',
    'v_name','j_name','c_name','source_contig_id','coding_frame','coding_aa_start','coding_aa_end',
    'coverage_kind','sequence_origin','source_nt_start','source_nt_end')


def load_inputs(bundle,geometry):
    """Strip every retained score/identity/correctness label before Evidence construction."""
    import polars as pl
    from vdjaqp.readout import _SCHEMA
    bundle=Path(bundle);candidate=bundle/'candidates.parquet'
    if not candidate.exists():candidate=bundle/f'{geometry}_candidate_identity.parquet'
    paths=[bundle/f'{geometry}_{enzyme}_psm_correctness.parquet' for enzyme in ('Trypsin','Chymotrypsin')]
    verify_bundle(bundle,[candidate,*paths])
    source_candidates=pl.read_parquet(candidate)
    # Source background-population provenance is independent of evaluation truth.
    columns=[*CANDIDATE_COLUMNS,*(['background_member'] if 'background_member' in source_candidates.columns else [])]
    candidates=source_candidates.select(columns)
    calls=[pl.read_parquet(p).select([*_SCHEMA,'enzyme']) for p in paths]
    peptides=pl.concat(calls,how='vertical_relaxed')
    if candidates['clone_id'].n_unique()!=candidates.height or peptides['peptide_id'].n_unique()!=peptides.height:
        raise ValueError('candidate/peptide identity is not unique')
    return candidates,peptides,[candidate,*paths]



def verify_bundle(bundle,paths):
    receipt=json.loads((Path(bundle)/'retained_table_manifest.json').read_text())
    artifacts={r['file']:r for r in receipt['artifacts']}
    for path in paths:
        name=str(Path(path).relative_to(bundle));record=artifacts.get(name)
        if record is None or path.stat().st_size!=record['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest()!=record['sha256']:
            raise ValueError('final scoring input differs from retained manifest: '+name)


def native_truth_pairs(metadata,peptides,truth_paths):
    """Spectrum × heavy anchored-junction class; all coisolated accession providers retained."""
    import polars as pl
    from vdjaqp.match import ms_collapse
    physical=metadata.filter(pl.col('locus')=='IGH')
    if physical['chain_id'].n_unique()!=physical.height:raise ValueError('physical chain IDs not unique')
    classes={r['chain_id']:s for r,s in zip(physical.to_dicts(),ms_collapse(physical['generator_junction_aa']))}
    universe=set(classes.values());all_ids=set(metadata['chain_id']);pairs=set();known=set()
    for enzyme,path in truth_paths.items():
        calls=peptides.filter(pl.col('enzyme')==enzyme)
        mapping={}
        for r in calls.select('native_spectrum_reference','spectrum_id').unique().to_dicts():
            if r['native_spectrum_reference'] in mapping:raise ValueError('native event maps to multiple observation keys')
            mapping[r['native_spectrum_reference']]=r['spectrum_id']
        path=Path(path)
        with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:events=json.load(f)
        seen=set()
        for event in events:
            key=event['native_spectrum_reference']
            if key in seen:raise ValueError('duplicate native event')
            seen.add(key)
            if key not in mapping:continue
            obs=mapping[key];known.add(obs)
            for contributor in event['contributors']:
                for accession in contributor['protein_accessions']:
                    if accession not in all_ids:raise ValueError('unknown native protein accession')
                    if accession in classes:pairs.add((obs,classes[accession]))
    return universe,known,pairs


def generative_n_geometry(metadata,generation):
    """Whole-codon TdT-only intervals from native generative records, never inferred NP."""
    from arda.refbuild.translate import translate
    raw={r['chain_id']:r for r in generation['chains']}
    if len(raw)!=len(generation['chains']) or set(raw)!=set(metadata['chain_id']):raise ValueError('generative/physical chain populations disagree')
    intervals={};unavailable={}
    for m in metadata.iter_rows(named=True):
        key=m['chain_id'];r=raw[key];a=r['airr_truth'];nt=r['variable_nt'].upper()
        if (r['protein_aa']!=m['protein_aa'] or r['variable_aa']!=m['variable_aa'] or
            translate(nt[:3*len(m['variable_aa'])])!=m['variable_aa'] or
            a['junction'].upper()!=m['generator_junction_nt'].upper() or a['junction_aa']!=m['generator_junction_aa'] or
            a['junction_start']!=m['generator_junction_start_nt'] or a['junction_end']!=m['generator_junction_end_nt'] or
            nt[a['junction_start']:a['junction_end']]!=a['junction'].upper()):
            raise ValueError('raw native generative record differs from physical metadata')
        if m['locus']!='IGH':continue
        bounds=[a.get(k) for k in ('v_sequence_end','d_sequence_start','d_sequence_end','j_sequence_start')]
        if any(type(x) is not int for x in bounds):unavailable[key]='missing_native_segment_geometry';continue
        ve,ds,de,js=bounds
        if not 0<=ve<=ds<=de<=js<=len(nt):raise ValueError('invalid native zero-based half-open segment geometry')
        parts=[];reason=None
        for number,(lo,hi) in enumerate(((ve,ds),(de,js)),1):
            name=f'np{number}';sequence=a.get(name)
            if not isinstance(sequence,str) or sequence.upper()!=nt[lo:hi] or a.get(name+'_length')!=hi-lo:
                raise ValueError('native NP slice/length differs from generative sequence')
            counts=[a.get(name+suffix) for suffix in ('_n_tdt','_n_p','_n_inserted')]
            if any(type(x) is not int or x<0 for x in counts):reason='missing_native_NP_origin_counts';break
            if sum(counts)!=hi-lo:raise ValueError('native NP origin counts do not sum to interval length')
            if counts[1] or counts[2]:reason='NP_contains_P_or_later_insertions_without_base_resolved_origin';break
            if not a['junction_start']<=lo<=hi<=a['junction_end']:raise ValueError('native N tract lies outside anchored junction')
            parts.append(((lo+2)//3,hi//3))
        if reason is None:intervals[key]=parts
        else:unavailable[key]=reason
    return intervals,unavailable


def native_junction_truth_pairs(metadata,peptides,truth_paths,generation,threads,placement_fn=None):
    """True contributor-peptide overlap with anchor-excluded physical heavy junctions."""
    import polars as pl
    from arda.refbuild.translate import translate
    from vdjaqp.match import ms_collapse
    if placement_fn is None:
        from vdjaqp.simulation_metrics import exact_peptide_placements
        placement_fn=exact_peptide_placements
    heavy=metadata.filter(pl.col('locus')=='IGH');geometry={};classes={};missing=set()
    n_geometry,n_unavailable=generative_n_geometry(metadata,generation)
    for r in heavy.iter_rows(named=True):
        key=r['chain_id'];junction=r['generator_junction_aa'];classes[key]=ms_collapse([junction])[0]
        start,end=r.get('generator_junction_start_nt'),r.get('generator_junction_end_nt')
        nt=r.get('generator_junction_nt')
        if start is None or end is None or nt is None:missing.add(key);continue
        if (type(start) is not int or type(end) is not int or start<0 or start%3 or end%3 or end<=start or
            end-start!=len(nt) or len(junction)<2 or junction[0]!='C' or junction[-1] not in 'FW' or
            translate(nt.upper())!=junction or r['protein_aa'][start//3:end//3]!=junction):
            raise ValueError('physical native junction coordinates/translation/protein slice disagree')
        geometry[key]=(start//3+1,end//3-1)
    records=[];known=set();queries={}
    all_accessions=set(metadata['chain_id'])
    for enzyme,path in truth_paths.items():
        calls=peptides.filter(pl.col('enzyme')==enzyme)
        mapping=dict(calls.select('native_spectrum_reference','spectrum_id').unique().iter_rows())
        with (gzip.open(path,'rt') if Path(path).suffix=='.gz' else Path(path).open()) as f:events=json.load(f)
        for event in events:
            if event['native_spectrum_reference'] not in mapping:continue
            obs=mapping[event['native_spectrum_reference']];known.add(obs)
            for contributor in event['contributors']:
                sequence=contributor['peptide_aa']
                for accession in contributor['protein_accessions']:
                    if accession not in all_accessions:raise ValueError('unknown native protein accession')
                    if accession in classes:
                        if len(sequence)>=6 and not set(sequence)-set('ACDEFGHIKLMNPQRSTVWY'):
                            queries.setdefault(sequence,str(len(queries)))
                        records.append((obs,accession,sequence))
    # One existing native whole-query batch, never one search per contributor/provider.
    q=pl.DataFrame({'peptide_id':list(queries.values()),'seq':list(queries)},schema={'peptide_id':pl.String,'seq':pl.String})
    hits=placement_fn(heavy.select('chain_id','protein_aa'),q,threads) if q.height else pl.DataFrame()
    proteins=dict(heavy.select('chain_id','protein_aa').iter_rows());sequences={v:k for k,v in queries.items()};positions={}
    for h in hits.iter_rows(named=True):
        accession=h['chain_id'];sequence=sequences[h['peptide_id']]
        if proteins[accession][h['start']:h['end']]==sequence:
            positions.setdefault((accession,sequence),set()).add((h['start'],h['end']))
    positive=set();unscorable=set();reasons={};n_positive=set();n_unscorable=set();n_reasons={}
    for obs,accession,sequence in records:
        spans=positions.get((accession,sequence),set())
        if accession not in geometry or not spans:
            unscorable.add(obs);reason='missing_physical_junction_geometry' if accession not in geometry else 'missing_exact_native_peptide_placement'
            reasons[reason]=reasons.get(reason,0)+1;n_unscorable.add(obs);n_reasons[reason]=n_reasons.get(reason,0)+1;continue
        left,right=geometry[accession]
        overlap={start<right and end>left for start,end in spans}
        if len(overlap)!=1:
            unscorable.add(obs);reasons['mixed_junction_framework_placements']=reasons.get('mixed_junction_framework_placements',0)+1
        elif True in overlap:positive.add((obs,classes[accession]))
        if accession not in n_geometry:
            n_unscorable.add(obs);reason=n_unavailable[accession];n_reasons[reason]=n_reasons.get(reason,0)+1
        else:
            n_overlap={any(start<right and end>left and left<right for left,right in n_geometry[accession]) for start,end in spans}
            if len(n_overlap)!=1:
                n_unscorable.add(obs);reason='mixed_N_nonN_placements';n_reasons[reason]=n_reasons.get(reason,0)+1
            elif True in n_overlap:n_positive.add((obs,classes[accession]))
    report=dict(native_observations=len(known),junction_scorable_observations=len(known-unscorable),
        junction_unscorable_observations=len(unscorable),unscorable_provider_record_reasons=reasons,
        heavy_proteins_with_junction_geometry=len(geometry),heavy_proteins_without_junction_geometry=len(missing),
        unique_native_heavy_contributor_peptides=len(queries),native_whole_query_batches=int(bool(queries)),
        strict_N_scorable_observations=len(known-n_unscorable),strict_N_unscorable_observations=len(n_unscorable),
        strict_N_unscorable_provider_record_reasons=n_reasons,heavy_proteins_with_verified_TdT_only_geometry=len(n_geometry),
        heavy_proteins_without_verified_TdT_only_geometry=len(n_unavailable),
        strict_generative_n_truth='at least one whole codon inside a verified native TdT-only tract; fixed labels across k; P/later-insertion tracts unavailable')
    return (known-unscorable,{(o,c) for o,c in positive if o not in unscorable},
        known-n_unscorable,{(o,c) for o,c in n_positive if o not in n_unscorable},report)


def attribution_metrics(predictions,universe,observations,positives,*,return_curve=False):
    """Exact pooled tied-score AP with implicit zero pairs; no all-present clone classification."""
    import numpy as np
    n=len(universe)*len(observations);p=len(positives)
    if any(o not in observations or c not in universe for o,c in positives):raise ValueError('truth outside universe')
    accum={}
    for o,c,score in predictions:
        if o not in observations or c not in universe or not np.isfinite(score) or score<0:
            raise ValueError('prediction outside evaluable universe')
        accum[o,c]=accum.get((o,c),0.)+float(score)
    nonzero={k:v for k,v in accum.items() if v>0}
    if any(v>1+1e-9 for v in nonzero.values()):raise ValueError('class mass exceeds observation mass')
    groups={}
    for key,score in nonzero.items():
        counts=groups.setdefault(score,[0,0]);counts[0]+=1;counts[1]+=key in positives
    tp=sum(k in positives for k in nonzero);fp=len(nonzero)-tp
    groups[0.]=[n-len(nonzero),p-tp]
    cumulative=positive_seen=0;ap=0.;points=[]
    for score,(count,pos) in sorted(groups.items(),reverse=True):
        if not count:continue
        cumulative+=count;positive_seen+=pos
        if p:ap+=pos/p*positive_seen/cumulative
        points.append(dict(threshold=score,true_positive=positive_seen,false_positive=cumulative-positive_seen,
            false_negative=p-positive_seen,precision=positive_seen/cumulative,recall=positive_seen/p if p else None,
            f1=2*positive_seen/(cumulative+p) if p else None,n_called=cumulative,positive_threshold=score>0))
    result=dict(n_observations=len(observations),n_junction_classes=len(universe),n_pairs=n,n_positive_pairs=p,
        nonzero_supported_pairs=len(nonzero),true_positive=tp,false_positive=fp,false_negative=p-tp,
        zero_score_positive_pairs=p-tp,
        precision=tp/(tp+fp) if tp+fp else 0.,recall=tp/p if p else None,
        f1=2*tp/(2*tp+fp+p-tp) if p else None,average_precision=ap if p else None,
        binary_call='strictly positive class explanation mass; not confirmed_nats>=20')
    eligible=[point for point in points if point['positive_threshold'] and point['f1'] is not None]
    # Points are descending by threshold, so equal F1 retains higher cutoff/fewer calls.
    best=max(eligible,key=lambda point:point['f1']) if eligible else dict(threshold=None,
        true_positive=0,false_positive=0,false_negative=p,precision=0.,recall=0. if p else None,f1=0. if p else None,n_called=0)
    result.update({'max_f1_'+key:best[key] for key in ('threshold','true_positive','false_positive','false_negative','precision','recall','f1','n_called')})
    rounded=rounded_operating_optimum(points,p)
    result.update({'rounded_max_f1_'+key:rounded[key] for key in ('threshold','true_positive','false_positive','false_negative','precision','recall','f1','n_called')})
    return (result,points) if return_curve else result


def rounded_operating_optimum(points,n_positive):
    """One sweep over two-significant-digit cutoff candidates; whole score ties stay intact."""
    from decimal import Decimal,ROUND_FLOOR
    cutoffs=set()
    for point in points:
        if not point['positive_threshold']:continue
        score=Decimal(str(point['threshold']));step=Decimal(1).scaleb(score.adjusted()-1)
        cutoffs.update((float(format(point['threshold'],'.2g')),float(score.quantize(step,rounding=ROUND_FLOOR))))
    cutoffs=sorted(cutoffs,reverse=True)
    empty=dict(true_positive=0,false_positive=0,false_negative=n_positive,precision=0.,
        recall=0. if n_positive else None,f1=0. if n_positive else None,n_called=0)
    best=dict(threshold=None,**empty);index=0;chosen=empty
    for cutoff in cutoffs:
        while index<len(points) and points[index]['positive_threshold'] and points[index]['threshold']>=cutoff:
            chosen=points[index];index+=1
        if chosen['f1'] is not None and (best['threshold'] is None or chosen['f1']>best['f1']):
            best=dict(threshold=cutoff,**{k:chosen[k] for k in empty})
    return best


def operating_point(points,threshold,n_positive):
    """Apply a fixed training threshold to complete held-out tie groups."""
    chosen=None
    if threshold is not None:
        for point in points:
            if point['positive_threshold'] and point['threshold']>=threshold:chosen=point
    chosen=chosen or dict(true_positive=0,false_positive=0,false_negative=n_positive,
        precision=0.,recall=0. if n_positive else None,f1=0. if n_positive else None,n_called=0)
    return {'operating_'+key:chosen[key] for key in ('true_positive','false_positive','false_negative','precision','recall','f1','n_called')}


def validate_operating_points(parameters):
    import math
    if (parameters.get('schema_version')!=1 or parameters.get('kind')!='synthetic_source_class_operating_points' or
        parameters.get('training_case')!='coding_only/M0-amplicon' or
        parameters.get('truth_definition')!='native_heavy_generative_N_full_codon_overlap'):
        raise ValueError('invalid operating-point model provenance/contract')
    for label in ('source_sha256','raw_input_sha256'):
        pins=parameters.get(label)
        if not isinstance(pins,dict) or not pins or any(not isinstance(v,str) or len(v)!=64 or set(v)-set('0123456789abcdef') for v in pins.values()):
            raise ValueError('operating-point model requires source/raw SHA256 provenance')
    points=parameters.get('operating_points',[]);seen=set()
    for point in points:
        k=point.get('minimum_n_run');threshold=point.get('threshold')
        if type(k) is not int or k in seen or (k!=0 and k<2):raise ValueError('invalid/duplicate operating-point N minimum')
        if threshold is not None and (type(threshold) not in (int,float) or not math.isfinite(threshold) or not 0<threshold<=1+1e-9 or threshold!=float(format(threshold,'.2g'))):raise ValueError('invalid rounded positive operating threshold')
        seen.add(k)
    if 0 not in seen or parameters.get('best_n_minimum') not in seen-{0}:raise ValueError('missing no-N or selected N operating point')
    return parameters


def learn_operating_points(rows,provenance):
    rows=[r for r in rows if r['aggregation']=='alignment_union' and r['view']!='all_observations_all_explanations' and
        r['truth_definition']=='native_heavy_generative_N_full_codon_overlap']
    points=[dict(minimum_n_run=r['minimum_n_run'],threshold=r['rounded_max_f1_threshold'],
        **{key:r['rounded_max_f1_'+key] for key in ('true_positive','false_positive','false_negative','precision','recall','f1','n_called')},
        average_precision=r['average_precision']) for r in rows]
    eligible=[r for r in points if r['minimum_n_run']>=2 and r['f1'] is not None]
    if not eligible:raise ValueError('M0 has no scorable strict-N training operating points')
    best=min(eligible,key=lambda r:(-r['f1'],-r['average_precision'],r['minimum_n_run']))
    return validate_operating_points(dict(schema_version=1,kind='synthetic_source_class_operating_points',
        training_case='coding_only/M0-amplicon',truth_definition='native_heavy_generative_N_full_codon_overlap',
        training_scope='fixed unexplained native-scorable observations and all physical junction classes',
        threshold_rule='class_score >= positive threshold; null means no calls; main 20-nat calls unchanged',
        cutoff_selection='maximum training F1 over two-significant-digit cutoffs; tie higher threshold/fewer calls',
        threshold_significant_digits=2,
        n_selection='maximum training F1, then AP, then smaller N minimum',
        source_sha256=provenance['source_sha256'],raw_input_sha256=provenance['raw_input_sha256'],
        operating_points=points,best_n_minimum=best['minimum_n_run']))


def rescue_evaluation(predictions,universe,observations,positives,mapped,*,return_curve=False):
    """One fixed unexplained universe for no-N and every N-run rescue view."""
    unexplained=observations-set(mapped)
    return attribution_metrics([(o,c,s) for o,c,s in predictions if o in unexplained],
        universe,unexplained,{(o,c) for o,c in positives if o in unexplained},return_curve=return_curve)


def annotate_candidates(candidates,threads,annotator=None):
    """Exactly one native cohort annotation; retain no Pgen/scenario probability bank."""
    import polars as pl
    if annotator is None:
        from vdjtools.model.junction import annotate_junctions
        annotator=annotate_junctions
    from vdjaqp.junction_support import observed_junction_classes
    observed=observed_junction_classes(candidates)
    eligible=candidates.filter(pl.col('clone_id').is_in(list(observed))&pl.col('v_name').is_not_null()&pl.col('j_name').is_not_null())
    fields=['cdr3_nt','v_end_nt','j_start_nt','d_start_nt','d_end_nt','np1','np2']
    if not eligible.height:return pl.DataFrame(schema={'clone_id':pl.String,**{k:pl.String if k in ('cdr3_nt','np1','np2') else pl.Int64 for k in fields}})
    native=annotator(eligible['cdr3_seq'].to_list(),eligible['v_name'].to_list(),eligible['j_name'].to_list(),species='human',threads=threads)
    if native.height!=eligible.height:raise ValueError('native annotation changed row population')
    native=native.with_columns(eligible['clone_id'])
    # Repaired sequences are not observed candidate sequences; retain unavailable geometry.
    same=native['cdr3_repaired']==eligible['cdr3_seq']
    return native.with_columns(*[pl.when(pl.Series(same)).then(pl.col(k)).otherwise(None).alias(k) for k in fields]).select('clone_id',*fields)


def j_scaffold_frames(germlines,markup,scaffolds):
    """Verify local J phase from exact public scaffold slices, never from an AA motif."""
    phases={allele:set() for allele in germlines}
    for row in markup.select('scaffold_id','j_call','j_sequence_start','vj_end','fwr4_start','fwr4_end','fwr4').iter_rows(named=True):
        try:start,end,anchor,stop=[int(row[k])-d for k,d in [('j_sequence_start',1),('vj_end',0),('fwr4_start',1),('fwr4_end',0)]]
        except (TypeError,ValueError):continue
        scaffold=scaffolds.get(row['scaffold_id'])
        if (scaffold is None or not 0<=start<=anchor<stop<=end<=len(scaffold) or
            not row['fwr4'] or scaffold[anchor:stop]!=row['fwr4']):continue
        for allele in (row['j_call'] or '').split(','):
            allele=allele.strip()
            if allele in phases and germlines[allele]==scaffold[start:end]:phases[allele].add((anchor-start)%3)
    qualified={a:next(iter(fs)) for a,fs in phases.items() if len(fs)==1}
    return qualified,{a:sorted(fs) for a,fs in phases.items() if len(fs)!=1}


def templated_competitors(peptides,threads):
    """One exact whole-query native batch against verified public V/J and forward-D ORFs."""
    import polars as pl
    from arda.germline import segment_germlines
    from arda.annotate.dmap import _d_germlines
    from arda.refbuild.translate import translate
    from arda.refbuild.imgt import read_fasta
    from arda.paths import vdj_dir
    from vdjaqp.match import PeptideMatcher
    reference=vdj_dir('human');rows=[]
    def add(name,region,frame,aa):
        offset=0
        for part in aa.split('*'):
            if part and set(part)<=set('ACDEFGHIKLMNPQRSTVWY') and len(part)>=6:
                rows.append(dict(template_id=f'{region}:{name}:f{frame}:aa{offset}',allele=name,region=region,
                    frame=frame,frame_aa_start=offset,seq=part,scope='public forward scaffold/reference; no inverse D'))
            offset+=len(part)+1
    for allele,nt in segment_germlines('human','v').items():
        if allele.startswith(('IGHV','IGKV','IGLV')):add(allele,'V',0,translate(nt))
    j={a:nt for a,nt in segment_germlines('human','j').items() if a.startswith(('IGHJ','IGKJ','IGLJ'))}
    frames,unavailable=j_scaffold_frames(j,pl.read_csv(reference/'markup.tsv',separator='\t',infer_schema_length=0),dict(read_fasta(reference/'alleles.fasta')))
    for allele,frame in frames.items():add(allele,'J',frame,translate(j[allele][frame:]))
    for allele,nt in _d_germlines('human').get('IGH',[]):
        for frame in range(3):add(allele,'D',frame,translate(nt[frame:]))
    templates=pl.DataFrame(rows)
    if not templates.height:raise ValueError('no qualified public templated references')
    hits=PeptideMatcher(templates['template_id'],templates['seq']).search(peptides['peptide_id'],peptides['seq'],max_subs=0,max_indels=0,threads=threads)
    hits=hits.rename({'clone_id':'template_id'}).join(templates.drop('seq'),on='template_id',validate='m:1')
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [reference/'alleles.fasta',reference/'markup.tsv',reference/'d_germlines.fasta']}
    return hits,templates,dict(raw_reference_sha256=hashes,j_without_unique_scaffold_phase=unavailable,
        template_counts_by_region_locus=templates.with_columns(pl.col('allele').str.slice(0,3).alias('locus')).group_by('region','locus').len().sort(['region','locus']).to_dicts(),
        scope='all returned public IGHV/IGKV/IGLV alleles in scaffold frame0; IGHJ/IGKJ/IGLJ with unanimous exact-scaffold coding phase; forward IGHD frames0,1,2; whole peptides only; not exhaustive biological germline diversity')


def compare_case(bundle,geometry,truth_paths,generation_path,output,threads,*,parameters,fit_operating,provenance):
    """Two fresh whole-cohort builds; all k views derive within each fresh build."""
    import gc,time
    import numpy as np
    import polars as pl
    from vdjaqp.pipeline import Evidence
    from vdjaqp.errormodel import ErrorModel,mass_prior,_E_DB_DEFAULT
    from vdjaqp.aggregate import ClonotypeScorer,spectrum_key
    from vdjaqp.junction_support import contextual_n_support,observed_junction_classes
    from vdjaqp.data import to_gene_level
    from vdjaqp.simulation_metrics import exact_peptide_placements,identity_metrics
    verify_bundle(Path(bundle),[Path(bundle)/'repertoire_metrics.parquet'])
    metadata=pl.read_parquet(Path(bundle)/'repertoire_metrics.parquet')
    if metadata.height!=2400 or metadata.filter(pl.col('locus')=='IGH').height!=1200:raise ValueError('complete physical population changed')
    output=Path(output);output.mkdir(parents=True,exist_ok=False);metrics=[];pair_receipt=None
    for mode in ('alignment_union','exclusive_spectrum'):
        mode_curves=[];curve_export=[]
        started=time.monotonic();dest=output/mode;dest.mkdir()
        candidates,peptides,inputs=load_inputs(bundle,geometry);candidates=to_gene_level(candidates)
        # Truth is loaded only by evaluation, never attached to candidate/call frames.
        universe,observations,positive=native_truth_pairs(metadata,peptides,truth_paths)
        with gzip.open(generation_path,'rt') as stream:generation=json.load(stream)
        junction_observations,junction_positive,n_observations,n_positive,junction_report=native_junction_truth_pairs(metadata,peptides,truth_paths,generation,threads)
        (dest/'junction_truth_scope.json').write_text(json.dumps(junction_report,indent=2)+'\n')
        model=ErrorModel(np.array([0.,1.]),np.array([.1,.01]),mass_prior(),n_substitutions=0,e_db=_E_DB_DEFAULT)
        ev=Evidence.from_frames(candidates,peptides,error_model=model)
        scores=ev.score(aggregation=mode);baseline_scores=ev.score(aggregation='alignment_union');scores.write_parquet(dest/'scores.parquet',compression='zstd')
        call_index=spectrum_key(peptides).select('peptide_id','spectrum_id','alc','alc_rank','n_reads')
        keyed=ev.alignments.join(call_index,on='peptide_id',validate='m:1')
        gamma=ClonotypeScorer._responsibility(keyed)
        keyed=keyed.with_columns(pl.Series('gamma',gamma))
        native=annotate_candidates(candidates,threads);native.write_parquet(dest/'inferred_geometry.parquet',compression='zstd')
        template_hits,templates,reference_report=templated_competitors(peptides,threads)
        template_hits.write_parquet(dest/'templated_alternatives.parquet',compression='zstd')
        templates.write_parquet(dest/'template_catalog.parquet',compression='zstd')
        (dest/'reference_scope.json').write_text(json.dumps(reference_report,indent=2)+'\n')
        credited=keyed.filter(pl.col('gamma')>0).select(pl.col('spectrum_id').alias('observation_id'),'clone_id').unique()
        placements=keyed.select(pl.col('spectrum_id').alias('observation_id'),'peptide_id','clone_id','start','end').with_columns(pl.lit(False).alias('is_germline'))
        support=contextual_n_support(candidates,peptides,placements,native,accepted_scores=baseline_scores.select('clone_id','confirmed_nats'),accepted_support=credited,min_consecutive=2)
        support=support.with_columns(pl.Series('gamma',gamma),pl.col('peptide_id').is_in(template_hits['peptide_id'].unique().implode()).alias('templated_alternative_present'))
        support.write_parquet(dest/'contextual_support.parquet',compression='zstd')
        receipt=dict(gamma_sha256=hashlib.sha256(keyed.select('clone_id','peptide_id','start','end','gamma').sort(['clone_id','peptide_id','start','end']).write_csv().encode()).hexdigest(),
            mapped_observations_sha256=hashlib.sha256(json.dumps(sorted(set(support.filter(pl.col('excluded_mapped_observation'))['observation_id']))).encode()).hexdigest(),
            evaluation_observations_sha256=hashlib.sha256(json.dumps(sorted(observations)).encode()).hexdigest(),
            junction_evaluation_sha256=hashlib.sha256(json.dumps([sorted(junction_observations),sorted(junction_positive)]).encode()).hexdigest(),
            strict_N_evaluation_sha256=hashlib.sha256(json.dumps([sorted(n_observations),sorted(n_positive)]).encode()).hexdigest())
        if pair_receipt is not None and receipt!=pair_receipt:raise ValueError('aggregation changed paired gamma, fixed-rule mask or observation universe')
        pair_receipt=receipt
        (dest/'paired_invariance.json').write_text(json.dumps(receipt,indent=2)+'\n')
        # Candidate junction class is an observed sequence key, not a truth feature.
        lookup=observed_junction_classes(candidates)
        represented_classes=set(lookup.values())
        rows=keyed.select('spectrum_id','clone_id','gamma').to_dicts()
        all_observations=set(call_index['spectrum_id'])
        mapped=set(support.filter(pl.col('excluded_mapped_observation'))['observation_id'])
        population_counts=dict(n_all_observations=len(all_observations),n_native_evaluable_observations=len(observations),
            n_explained_native_observations=len(observations&mapped),n_unexplained_native_observations=len(observations-mapped),
            n_unscorable_observations=len(all_observations-observations),n_mapped_all_observations=len(mapped))
        baseline=[(r['spectrum_id'],lookup[r['clone_id']],r['gamma']) for r in rows
                  if r['spectrum_id'] in observations and lookup.get(r['clone_id']) in universe]
        prediction_views=[(0,'all_observations_all_explanations',baseline),(0,'unexplained_rescue_no_n_filter',baseline)]
        maxrun=int(support['longest_n_run'].max() or 0)
        grid=set(range(2,max(3,maxrun+1)+1))
        if parameters is not None:grid.update(p['minimum_n_run'] for p in parameters['operating_points'] if p['minimum_n_run']>=2)
        for k in sorted(grid):
            selected=support.filter((pl.col('longest_n_run')>=k)&~pl.col('excluded_mapped_observation'))
            pred=[(r['observation_id'],lookup[r['clone_id']],r['gamma']) for r in selected.select('observation_id','clone_id','gamma').to_dicts()
                  if r['observation_id'] in observations and lookup.get(r['clone_id']) in universe]
            prediction_views.append((k,'unexplained_rescue',pred))
        for truth_definition,evaluable,positives in [
            ('any_native_heavy_provider',observations,positive),
            ('native_heavy_junction_interior_overlap',junction_observations,junction_positive),
            ('native_heavy_generative_N_full_codon_overlap',n_observations,n_positive)]:
            counts=dict(n_all_observations=len(all_observations),n_native_evaluable_observations=len(evaluable),
                n_explained_native_observations=len(evaluable&mapped),n_unexplained_native_observations=len(evaluable-mapped),
                n_unscorable_observations=len(all_observations-evaluable),n_mapped_all_observations=len(mapped))
            for k,view,pred in prediction_views:
                rescue=view!='all_observations_all_explanations'
                result,curve=(rescue_evaluation(pred,universe,evaluable,positives,mapped,return_curve=True) if rescue else
                    attribution_metrics([(o,c,v) for o,c,v in pred if o in evaluable],universe,evaluable,positives,return_curve=True))
                unrepresented=sum(c not in represented_classes and (not rescue or o not in mapped) for o,c in positives)
                metrics.append(dict(aggregation=mode,minimum_n_run=k,view=view,truth_definition=truth_definition,
                    positive_pairs_without_observed_candidate=unrepresented,**counts,**result))
                mode_curves.append((len(metrics)-1,curve))
                curve_export.extend(dict(minimum_n_run=k,view=view,truth_definition=truth_definition,**point) for point in curve)
        if fit_operating and parameters is None:parameters=learn_operating_points(metrics,provenance)
        cutoffs={point['minimum_n_run']:point['threshold'] for point in parameters['operating_points']}
        for index,curve in mode_curves:
            row=metrics[index];k=row['minimum_n_run'];available=k in cutoffs
            row.update(operating_point_available=available,operating_threshold=cutoffs.get(k),
                selected_training_N=k==parameters['best_n_minimum'],operating_point_source='original M0 strict-N training only')
            if available:row.update(operating_point(curve,cutoffs[k],row['n_positive_pairs']))
            # A held-out maximum is not promoted or exported as an operating point.
            if not (fit_operating and row['truth_definition']=='native_heavy_generative_N_full_codon_overlap' and row['view']!='all_observations_all_explanations'):
                for key in list(row):
                    if key.startswith(('max_f1_','rounded_max_f1_')):del row[key]
        pl.DataFrame(curve_export,strict=False).write_csv(dest/'source_class_precision_recall.tsv',separator='\t')
        unrepresented_all=sum(c not in represented_classes for o,c in positive)
        unrepresented_rescue=sum(o not in mapped and c not in represented_classes for o,c in positive)
        heavy=metadata.filter(pl.col('locus')=='IGH').select('chain_id','protein_aa')
        query=candidates.select(pl.col('clone_id').alias('peptide_id'),'seq')
        assoc=exact_peptide_placements(heavy,query,threads).rename({'peptide_id':'clone_id'})
        ranks=scores.sort(['score','clone_id'],descending=[True,False]).with_row_index('rank',offset=1)
        population,called,identity=identity_metrics(metadata,candidates,ranks,assoc)
        for name,frame in [('identity_population',population),('candidate_identity',called),('identity_metrics',identity)]:frame.write_parquet(dest/(name+'.parquet'),compression='zstd')
        report=dict(aggregation=mode,n_candidates=candidates.height,n_peptides=peptides.height,n_physical_heavy=heavy.height,
            n_placements=keyed.height,n_native_evaluable_observations=len(observations),n_positive_pairs=len(positive),n_junction_classes=len(universe),
            unresolved_junction_placements=sum(lookup.get(r['clone_id']) not in universe for r in rows),
            unresolved_junction_mass=sum(r['gamma'] for r in rows if lookup.get(r['clone_id']) not in universe),
            seconds=time.monotonic()-started,minimum_grid=sorted(grid),
            source_input_sha256={str(p):hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in inputs},
            truth_used_for_scoring=False,standalone_germline_role='exact reference annotation only; not extra gamma component',
            candidate_NP_origin='model-inferred non-template tract; exported boundaries do not prove TdT-only origin',
            source_ap_aggregation_invariant_by_definition=True,n_geometry_missing=support.filter(pl.col('geometry_status')=='missing_geometry').height,
            n_candidates_with_complete_concordant_junction=len(lookup),
            n_candidates_with_available_np_geometry=native.filter(pl.col('cdr3_nt').is_not_null()&pl.col('v_end_nt').is_not_null()&pl.col('j_start_nt').is_not_null()&pl.col('d_start_nt').is_not_null()&pl.col('d_end_nt').is_not_null()).height,
            observation_population=population_counts,
            all_observation_positive_pairs_without_observed_candidate=unrepresented_all,
            rescue_positive_pairs_without_observed_candidate=unrepresented_rescue,
            n_global_peptide_call_candidate_record_pairs=peptides.height*candidates.height,
            global_collision_family='all peptide/repertoire comparisons; no empirical collision/null computation here')
        (dest/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        del ev,keyed,support;gc.collect()
    pl.DataFrame(metrics,strict=False).write_csv(output/'attribution_metrics.tsv',separator='\t')
    return metrics,parameters



def file_hashes(paths):
    return {str(Path(p).resolve()):hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}


def verify_pins(pins):
    if file_hashes(pins)!=pins:raise ValueError('pinned source/input/reference files changed during comparison')


def resolve_generation_truth(cases,leader_root):
    """Reuse existing raw generation files only after exact original producer-hash matching."""
    leader_root=Path(leader_root);available={}
    paths=[*leader_root.glob('M*-rnaseq/truth.json.gz'),*(leader_root/'original_truth').rglob('*.json.gz')]
    for path in sorted(paths):
        with gzip.open(path,'rb') as stream:sha=hashlib.file_digest(stream,'sha256').hexdigest()
        available.setdefault(sha,path)
    selected={};receipt=[]
    for model,name,bundle,_,_ in cases:
        producer=Path(bundle)/'metric_integration.json'
        sha=json.loads(producer.read_text())['artifact_sha256']['truth.json']
        if sha not in available:raise ValueError('original generative raw truth missing: '+model+'/'+name)
        path=available[sha];selected[model,name]=path
        receipt.append(dict(model=model,case=name,path=str(path.resolve()),raw_sha256=sha,
            gzip_sha256=file_hashes([path])[str(path.resolve())],producer_receipt=str(producer.resolve())))
    return selected,receipt


def selected_case_indices(case_index=None,case_indices=None):
    """Explicit fresh cohort selection; never infer selection from previous outputs."""
    if case_index is not None and case_indices is not None:raise ValueError('case-index and case-indices are mutually exclusive')
    selected=list(range(18)) if case_index is None and case_indices is None else ([case_index] if case_index is not None else list(case_indices))
    if not selected or len(set(selected))!=len(selected) or any(type(i) is not int or not 0<=i<18 for i in selected):
        raise ValueError('selected indices must be unique integers in 0..17')
    return selected


def main():
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original-root',type=Path,required=True);p.add_argument('--leader-root',type=Path,required=True)
    p.add_argument('--truth-manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--operating-points',type=Path,help='fixed M0-learned model parameters; required when M0 is not selected')
    p.add_argument('--cpus',type=int,default=16)
    selection=p.add_mutually_exclusive_group()
    selection.add_argument('--case-index',type=int)
    selection.add_argument('--case-indices',type=int,nargs='+')
    a=p.parse_args();threads=resolve_cpu_budget(a.cpus)
    try:selected_indices=selected_case_indices(a.case_index,a.case_indices)
    except ValueError as error:p.error(str(error))
    if (0 in selected_indices)==(a.operating_points is not None):p.error('fit M0 without operating-points, or supply fixed operating-points for a cohort excluding M0')
    parameters=validate_operating_points(json.loads(a.operating_points.read_text())) if a.operating_points else None
    import polars as pl
    truth=json.loads(a.truth_manifest.read_text())
    if truth['status']!='complete' or truth['retrieved']!=36:p.error('all thirty-six verified truth files required')
    rows=truth['files'];cases=[]
    for model,folder,geometries in [('coding_only',a.original_root,('amplicon','rnaseq')),('leader_imputed_v1',a.leader_root,('rnaseq',))]:
        for number in range(6):
            for geometry in geometries:
                name=f'M{number}-{geometry}';bundle=folder/name
                files={r['enzyme']:r for r in rows if r['model']==model and r['case']==name}
                if set(files)!={'Trypsin','Chymotrypsin'}:p.error('incomplete per-case truth')
                paths={}
                for enzyme,r in files.items():
                    path=Path(r['local'])
                    if hashlib.sha256(path.read_bytes()).hexdigest()!=r['compressed_sha256']:p.error('truth compressed hash mismatch')
                    paths[enzyme]=path
                cases.append((model,name,bundle,geometry,paths))
    if len(cases)!=18:p.error('invalid complete inventory')
    if a.output.exists():p.error('fresh output required; no resume')
    generation_paths,generation_receipt=resolve_generation_truth(cases,a.leader_root)
    import arda,vdjtools,seqtree
    from arda.paths import vdj_dir
    repo=Path(__file__).resolve().parents[1];reference=vdj_dir('human')
    source_paths=[*sorted((repo/'vdjaqp').glob('*.py')),Path(__file__).resolve()]
    arda_root=Path(arda.__file__).parent;vdj_root=Path(vdjtools.__file__).parent
    source_paths += [arda_root/'germline.py',arda_root/'cdr3fix.py',arda_root/'annotate'/'dmap.py',
        arda_root/'refbuild'/'build.py',arda_root/'refbuild'/'combinations.py',
        *sorted((vdj_root/'model').glob('*.py')),Path(seqtree.__file__)]
    raw_paths=[a.truth_manifest,*generation_paths.values(),*[p for _,_,_,_,paths in cases for p in paths.values()]]
    if a.operating_points:raw_paths.append(a.operating_points)
    for _,_,bundle,geometry,_ in cases:
        candidate=bundle/'candidates.parquet'
        if not candidate.exists():candidate=bundle/f'{geometry}_candidate_identity.parquet'
        raw_paths.extend([bundle/'retained_table_manifest.json',bundle/'metric_integration.json',bundle/'repertoire_metrics.parquet',candidate,
            *[bundle/f'{geometry}_{enzyme}_psm_correctness.parquet' for enzyme in ('Trypsin','Chymotrypsin')]])
    reference_paths=[reference/n for n in ('alleles.fasta','markup.tsv','d_germlines.fasta','cdr3_anchors.tsv')]
    for source in ('olga','arda'):
        for locus in ('IGH','IGK','IGL'):
            key=f'human_{locus}' if source=='arda' else locus
            reference_paths.extend(p for p in (vdj_root/'model'/'_bundled'/source/key).rglob('*') if p.is_file())
    provenance=dict(status='running',cpu_budget=threads,source_sha256=file_hashes(source_paths),
        raw_input_sha256=file_hashes(raw_paths),public_reference_sha256=file_hashes(reference_paths),
        native_generative_truth=generation_receipt,
        primary_truth='native heavy provider peptide overlaps at least one whole codon in verified generative TdT-only N',
        requested_cases=len(selected_indices),qualification_index=a.case_index,selected_case_indices=selected_indices,
        operating_point_input=str(a.operating_points.resolve()) if a.operating_points else None,
        truth_used_for_scoring=False,native_geometry='inferred, not measured; fresh batch per Evidence configuration')
    a.output.mkdir(parents=True);all_metrics=[]
    def save(): (a.output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    save()
    try:
        for index,(model,name,bundle,geometry,paths) in enumerate(cases):
            if index not in selected_indices:continue
            results,parameters=compare_case(bundle,geometry,paths,generation_paths[model,name],a.output/model/name,threads,
                parameters=parameters,fit_operating=index==0,provenance=provenance)
            if index==0:(a.output/'source_class_operating_points.json').write_text(json.dumps(parameters,indent=2,allow_nan=False)+'\n')
            all_metrics.extend(dict(model=model,case=name,**r) for r in results)
            pl.DataFrame(all_metrics,strict=False).write_csv(a.output/'metrics.tsv',separator='\t')
        # Tuning selection is confined to original M0 amplicon and retains the entire grid.
        training=pl.DataFrame(all_metrics,strict=False).filter((pl.col('model')=='coding_only')&(pl.col('case')=='M0-amplicon')&(pl.col('view')!='all_observations_all_explanations')&(pl.col('truth_definition')=='native_heavy_generative_N_full_codon_overlap'))
        training.write_csv(a.output/'training_parameter_grid.tsv',separator='\t')
        pareto=[]
        for row in training.to_dicts():
            if row['average_precision'] is None or row['operating_f1'] is None:continue
            dominated=any(other['average_precision'] is not None and other['operating_f1'] is not None and
                other['average_precision']>=row['average_precision'] and other['operating_f1']>=row['operating_f1'] and
                (other['average_precision']>row['average_precision'] or other['operating_f1']>row['operating_f1'])
                for other in training.to_dicts())
            if not dominated:pareto.append(row)
        pl.DataFrame(pareto,schema=training.schema).write_csv(a.output/'training_pareto.tsv',separator='\t')
        for label in ('source_sha256','raw_input_sha256','public_reference_sha256'):verify_pins(provenance[label])
        (a.output/'execution.json').write_text(json.dumps(dict(status='complete',cases=len({(r['model'],r['case']) for r in all_metrics}),
            qualification_index=a.case_index,fresh_evidence_builds=2*len({(r['model'],r['case']) for r in all_metrics}),
            selected_case_indices=selected_indices,
            tuning_case='coding_only/M0-amplicon strict-generative-N only',other_cases_untuned=True,source_ap_does_not_test_aggregation=True),indent=2)+'\n')
        provenance['status']='complete'
    except BaseException as error:
        provenance.update(status='failed',error=repr(error));raise
    finally:
        provenance['artifact_sha256']=file_hashes(p for p in a.output.rglob('*') if p.is_file() and p.name!='provenance.json')
        save()

if __name__=='__main__':main()
