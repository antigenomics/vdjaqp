"""Explicit source-region, acquisition and pooled evaluation stages for a matched mock.

No generation/resume; truth is evaluation-only. One physical solution trypsin,
one heavy-band trypsin and one shared heavy chymotrypsin per geometry.
"""
import argparse
from datetime import datetime,timezone
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.export_mock_observation_metrics import digest,raw_truth_digest,verify_observations,bounded_coverage

BLOCKS=('solution_trypsin','heavy_band_trypsin','shared_chymotrypsin')


def acquisition_scope(block):
    if block not in BLOCKS:raise ValueError('unknown acquisition block')
    return ('Chymotrypsin' if block==BLOCKS[2] else 'Trypsin',block!=BLOCKS[0])


def read_fasta(path):
    records={};key=None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            key=line[1:].split()[0]
            if key in records:raise ValueError('duplicate FASTA identity')
            records[key]=''
        elif line:
            if key is None:raise ValueError('sequence before FASTA identity')
            records[key]+=line.strip()
    if any(not seq or set(seq)-set('ACDEFGHIKLMNPQRSTVWY') for seq in records.values()):raise ValueError('noncanonical FASTA protein')
    return records


def verify_manifest(path):
    record=json.loads(path.read_text())
    if record['status']!='complete':raise ValueError('incomplete input manifest')
    for name,h in record['artifact_sha256'].items():
        f=(path.parent/name).resolve()
        if not f.is_relative_to(path.parent.resolve()) or digest(f)!=h:raise ValueError('input artifact hash mismatch')
    return record


def packet_frames(path):
    import polars as pl
    record=verify_manifest(path);base=path.parent
    for name,h in record['source_input_sha256'].items():
        if digest(Path(name))!=h:raise ValueError('qualified source pin changed')
    sources=json.loads((base/'reported_sources.json').read_text())['records']
    proteins=read_fasta(base/'physical_paired_reference_model.fasta');decisions={d['reference']:d for d in record['decisions']}
    panel=pl.read_parquet(base/'supplied_heavy_query_panel.parquet');rows=[]
    for r in sources:
        for side in ('H','L'):
            chain=r['reference']+'_'+side;variable=r['heavy_aa' if side=='H' else 'light_aa'];protein=proteins[chain]
            if not protein.startswith(variable):raise ValueError('reported variable/physical model mismatch')
            query=panel.filter(pl.col('clone_id')==r['clone_id']).row(0,named=True)
            rows.append(dict(chain_id=chain,locus='IGH' if side=='H' else decisions[r['reference']]['chosen_model_locus'],variable_aa=variable,protein_aa=protein,variable_length=len(variable),reference=r['reference'],query_clone_id=r['clone_id'] if side=='H' else None,known_junction=query['cdr3_seq'] if side=='H' else None,constant_gene='IGHG1' if side=='H' else decisions[r['reference']]['constant_gene'],provided_v_call=query['v_name'] if side=='H' else r['reported_light_v'],provided_j_call=query['j_name'] if side=='H' else r['reported_light_j']))
    metadata=pl.DataFrame(rows)
    if panel.height!=12 or metadata.height!=24 or metadata['chain_id'].n_unique()!=24:raise ValueError('expected twelve paired references')
    return record,panel,metadata


def aa_regions(metadata,annotations):
    """ARDA AA coordinates are 1-based inclusive; never invent NT/genotype truth."""
    import polars as pl
    by_id={a['sequence_id']:a for a in annotations}
    if len(by_id)!=len(annotations) or set(by_id)!=set(metadata['chain_id']):raise ValueError('AA annotation population mismatch')
    rows=[]
    for m in metadata.iter_rows(named=True):
        a=by_id[m['chain_id']];var=m['variable_aa'];n=len(var)
        if a['sequence']!=var:raise ValueError('annotation sequence changed')
        junction=m.get('known_junction') or a.get('junction_aa')
        if junction and (var.count(junction)!=1 or junction[0]!='C' or junction[-1] not in 'FW'):raise ValueError('ambiguous/nonanchored source junction')
        junction_start=var.index(junction) if junction else None
        def add(region,start,end,source,gene=None,status=None):
            if start is not None and not 0<=start<end<=len(m['protein_aa']):raise ValueError('invalid AA region')
            rows.append(dict(chain_id=m['chain_id'],locus=m['locus'],region=region,start=start,end=end,status=status or ('available' if start is not None else 'unavailable'),source=source,gene_call=gene))
        for native,label in [('fwr1','FR1'),('cdr1','CDR1'),('fwr2','FR2'),('cdr2','CDR2'),('fwr3','FR3'),('cdr3','CDR3'),('fwr4','FR4')]:
            lo,hi=a.get(native+'_start'),a.get(native+'_end')
            lo,hi=(int(lo)-1,int(hi)) if lo not in (None,'') and hi not in (None,'') else (None,None)
            if lo is not None and (hi>n or a.get(native+'_aa') not in (None,'',var[lo:hi])):raise ValueError('AA coordinate/slice mismatch')
            conflict=False
            if junction and lo is not None:
                stop=junction_start+len(junction)
                conflict=((label=='FR3' and hi!=junction_start+1) or (label=='CDR3' and (lo!=junction_start+1 or hi!=stop-1)) or (label=='FR4' and lo!=stop-1))
            add(label,None if conflict else lo,None if conflict else hi,'native_AA_inferred',status='unavailable_junction_boundary_conflict' if conflict else None)
        if junction:add('JUNCTION',junction_start,junction_start+len(junction),'supplied_or_native_AA_anchored')
        else:add('JUNCTION',None,None,'source_junction_unavailable')
        add('VARIABLE',0,n,'reported_variable');add('PROTEIN',0,len(m['protein_aa']),'declared_physical_model')
        add('C',n,len(m['protein_aa']),'public_constant_model',m.get('constant_gene'))
        for segment,lo,hi in [('V',a.get('v_sequence_start'),a.get('v_sequence_end')),('J',a.get('j_sequence_start'),a.get('fwr4_end'))]:
            start,end=(int(lo)-1,int(hi)) if lo not in (None,'') and hi not in (None,'') else (None,None)
            if end is not None and end>n:raise ValueError('inferred segment outside variable')
            add(segment,start,end,'native_AA_inferred_segment',a.get(segment.lower()+'_call'))
        for region in ('D','NP1','NP2'):add(region,None,None,'unavailable_from_reported_AA')
    return pl.DataFrame(rows).with_columns(pl.col('start','end').cast(pl.Int64))


def class_benchmark(universe,full_pairs,scores):
    """True full-variable classes; only unambiguous full candidate mappings credit scores."""
    import polars as pl
    from vdjaqp.bench import classification_metrics
    if scores['clone_id'].null_count() or scores['clone_id'].n_unique()!=scores.height or scores['score'].null_count() or scores['confirmed_nats'].null_count() or not scores['score'].is_finite().all() or not scores['confirmed_nats'].is_finite().all():raise ValueError('invalid candidate score population')
    if full_pairs.join(scores.select('clone_id'),on='clone_id',how='anti').height or full_pairs.join(universe.select('variable_mass_class_id').unique(),on='variable_mass_class_id',how='anti').height:raise ValueError('unknown evaluation association identity')
    u=universe.group_by('variable_mass_class_id').agg(pl.col('physical_positive').any())
    pairs=full_pairs.select('clone_id','variable_mass_class_id').unique()
    mapping=scores.join(pairs.group_by('clone_id').agg(pl.len().alias('full_class_count'),pl.col('variable_mass_class_id').first()),on='clone_id',how='left').with_columns(pl.col('full_class_count').fill_null(0))
    mapping=mapping.with_columns(pl.when(pl.col('full_class_count')==1).then(pl.col('variable_mass_class_id')).otherwise(None).alias('variable_mass_class_id'),pl.when(pl.col('full_class_count')==1).then(pl.lit('unique_full_class')).when(pl.col('full_class_count')==0).then(pl.lit('no_full_class')).otherwise(pl.lit('ambiguous_full_classes')).alias('assignment_status'))
    observed=(mapping.filter(pl.col('full_class_count')==1).group_by('variable_mass_class_id').agg(pl.col('score').max(),pl.col('confirmed_nats').max()))
    minimum=observed['score'].min();sentinel=minimum-1 if minimum is not None else -1.
    if minimum is not None and not sentinel<minimum:raise ValueError('cannot represent below-observed evaluation sentinel')
    classes=(u.join(observed,on='variable_mass_class_id',how='left').join(pairs.select('variable_mass_class_id').unique().with_columns(pl.lit(True).alias('database_observable')),on='variable_mass_class_id',how='left')
             .with_columns(pl.col('score').fill_null(sentinel),pl.col('confirmed_nats').fill_null(0.),pl.col('database_observable').fill_null(False)).rename({'variable_mass_class_id':'clone_id'}).sort('clone_id'))
    positive=set(classes.filter(pl.col('physical_positive'))['clone_id']);metrics={'end_to_end_full_variable_IL_classes':classification_metrics(classes,positive)}
    conditional=classes.filter(pl.col('database_observable'));cp=positive&set(conditional['clone_id'])
    metrics['database_observable_full_variable_IL_classes']=classification_metrics(conditional,cp) if cp else dict(status='no_observable_positive_classes',n_candidates=conditional.height,n_positives=0)
    called=mapping.filter(pl.col('confirmed_nats')>=20);tp=called.filter((pl.col('full_class_count')==1)&pl.col('variable_mass_class_id').is_in(positive)).height
    unresolved=called.filter(pl.col('full_class_count')!=1).height
    metrics['all_shortlist_candidate_records']=dict(called_records=called.height,definitive_positive_records=tp,unresolved_called_records=unresolved,definitive_negative_records=called.height-tp-unresolved,definitive_precision_all_calls=tp/called.height if called.height else 0.)
    metrics['evaluation_only_missing_class_score']=sentinel
    return classes,mapping.sort('clone_id'),metrics


def build_database(entries,light_aux,constants,path):
    names=set(entries['entry_id'])
    if len(names)!=entries.height or names&set(light_aux) or (names|set(light_aux))&set(constants):raise ValueError('search entry identity collision')
    with path.open('x') as f:
        for name,seq in [*entries.select('entry_id','seq').iter_rows(),*light_aux.items(),*constants.items()]:
            if not seq or set(seq)-set('ACDEFGHIKLMNPQRSTVWY'):raise ValueError('invalid database sequence')
            f.write('>'+name+'\n'+seq+'\n')


def retire_native(native,*,validated):
    if not validated:raise ValueError('native exports not validated')
    for name in ('psms.idXML','peptides.tsv','ms2_truth.json'):
        if not (native/name).is_file():raise ValueError('missing final native readout')
    ledger=[]
    for name in ('spectra.mzML','truth.featureXML'):
        path=native/name
        if not path.resolve().is_relative_to(native.resolve()) or not path.is_file():raise ValueError('unsafe/missing owned native artifact')
        record=dict(path=name,raw_sha256=digest(path),bytes=path.stat().st_size,validated_at=datetime.now(timezone.utc).isoformat());path.unlink();record['status']='deleted_after_validation';ledger.append(record)
    return ledger


def write_frames(output,frames):
    for name,frame in frames.items():frame.write_parquet(output/(name+'.parquet'),compression='zstd')



def input_pins(paths):
    repo=Path(__file__).resolve().parents[1]
    code=[Path(__file__),*sorted((repo/'vdjaqp').glob('*.py')),*[repo/'scripts'/n for n in ['smoke_proteomics.py','export_mock_observation_metrics.py','smoke_metrics.py']]]
    return {str(p.resolve()):digest(p) for p in [*paths,*code]}


def complete_manifest(output,record,inputs):
    if any(digest(Path(f))!=h for f,h in inputs.items()):raise ValueError('input changed during stage')
    record.update(status='complete',input_sha256=inputs,source_sha256={str(f.resolve()):digest(f) for f in [Path(__file__),*sorted((Path(__file__).resolve().parents[1]/'vdjaqp').glob('*.py'))]},artifact_sha256={str(f.relative_to(output)):digest(f) for f in output.rglob('*') if f.is_file()})
    (output/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')


def qualify_regions(args):
    from arda.annotate.contig import reannotate_contigs
    _,_,metadata=packet_frames(args.packet_manifest)
    inputs=input_pins([args.packet_manifest])
    # One native batch of the actual reported VAA, never reverse-translated NT.
    annotations=reannotate_contigs(list(metadata.select('chain_id','variable_aa').iter_rows()),seqtype='aa',threads=args.cpus,map_d=False)
    regions=aa_regions(metadata,annotations)
    import polars as pl
    calls=pl.DataFrame([dict(chain_id=a['sequence_id'],v_call=a.get('v_call') or None,j_call=a.get('j_call') or None) for a in annotations])
    metadata=metadata.join(calls,on='chain_id',how='left',validate='1:1').with_columns(pl.col('constant_gene').alias('constant_nt_gene'),pl.lit(None,dtype=pl.String).alias('canonical_v_gene'),pl.lit('native_AA_inferred_not_canonical_mapped').alias('canonical_v_match_status'),pl.lit(None,dtype=pl.Int64).alias('canonical_v_gene_count'))
    args.output.mkdir(parents=True)
    (args.output/'native_AA_annotations.json').write_text(json.dumps(annotations,indent=2)+'\n')
    write_frames(args.output,{'physical_metadata':metadata,'physical_regions':regions})
    complete_manifest(args.output,dict(stage='source_AA_regions',coordinate_units='ARDA AA one-based inclusive converted to zero-based half-open',chains=metadata.height,annotation_batches=1,D_N_genealogy='unavailable; no generated genotype',packet_manifest_sha256=digest(args.packet_manifest)),inputs)


def inference_inputs(observed,packet_manifest):
    import polars as pl
    from vdjaqp.reconstruction import append_known_references
    record,panel,metadata=packet_frames(packet_manifest)
    candidates=pl.read_parquet(observed/'candidates.parquet');entries=pl.read_parquet(observed/'search_entries.parquet')
    result=append_known_references(candidates,entries,panel,expected_count=12,source_format='retained supplied heavy-query packet; original fields preserved')
    return result,metadata


def acquire(args):
    import polars as pl
    from vdjaqp.readout import read_idxml
    from vdjaqp.simulation_metrics import psm_correctness
    verify_manifest(args.reference_regions)
    if json.loads(args.reference_regions.read_text())['packet_manifest_sha256']!=digest(args.packet_manifest):raise ValueError('physical region/packet mismatch')
    _,_,observed=verify_observations(args.observed_receipt,args.geometry)
    inference,physical=inference_inputs(observed,args.packet_manifest)
    enzyme,heavy_only=acquisition_scope(args.block)
    selected=physical.filter(pl.col('locus')=='IGH') if heavy_only else physical
    args.output.mkdir(parents=True)
    database=args.output/'observed_plus_supplied_search.fasta'
    build_database(inference['search_entries'],read_fasta(args.packet_manifest.parent/'supplied_light_auxiliary_search.fasta'),read_fasta(args.packet_manifest.parent/'public_constant_search.fasta'),database)
    physical_path=args.output/'physical_model_input.fasta'
    with physical_path.open('x') as f:
        for name,seq in selected.select('chain_id','protein_aa').iter_rows():f.write('>'+name+' [# intensity=10000 #]\n'+seq+'\n')
    write_frames(args.output,inference)
    inputs=input_pins([args.packet_manifest,args.reference_regions,args.observed_receipt,args.mssimulator])
    native=args.output/'native';run_id=args.geometry+':'+args.block
    command=[sys.executable,str(Path(__file__).with_name('smoke_proteomics.py')),str(physical_path.resolve()),str(database.resolve()),str(native.resolve()),'--enzyme',enzyme,'--threads',str(args.cpus),'--gradient-seconds','1200','--ms2-per-bin','5','--min-detectability','0.5','--dynamic-exclusion-seconds','0','--keep-spectra','--mssimulator',str(args.mssimulator.resolve()),'--stage-timeout',str(args.stage_timeout)]
    with (args.output/'native.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=2*args.stage_timeout+60)
    calls=read_idxml(native/'psms.idXML',run_id=run_id)
    native_truth=json.loads((native/'ms2_truth.json').read_text())
    if any(set(c['protein_accessions'])-set(selected['chain_id']) for e in native_truth for c in e['contributors']):raise ValueError('native truth names unassayed protein')
    per,accuracy=psm_correctness(calls,native_truth)
    if per.filter(pl.col('status')!='evaluated').height:raise ValueError('unlinked native PSM assignments')
    write_frames(args.output,{'raw_calls':calls,'psm_correctness':per,'psm_accuracy':accuracy})
    # Producer summary remains byte-intact. Retirement records actual owned paths.
    before=sum(f.stat().st_size for f in args.output.rglob('*') if f.is_file())
    ledger=retire_native(native,validated=True)
    after=sum(f.stat().st_size for f in args.output.rglob('*') if f.is_file())
    retirement=dict(files=ledger,task_logical_bytes_before=before,task_logical_bytes_after=after,whole_project_disk='outer process-tree profiler receipt; task bytes are not project allocation')
    (args.output/'retirement.json').write_text(json.dumps(retirement,indent=2)+'\n')
    # Requested raw native finals are losslessly retained, never execution features.
    retention=[]
    for name in ('psms.idXML','truth.idXML','peptides.tsv','ms2_truth.json'):
        src=native/name;packed=src.with_suffix(src.suffix+'.gz');raw_hash=digest(src);raw_bytes=src.stat().st_size
        with src.open('rb') as inp,packed.open('xb') as raw:
            with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as dest:shutil.copyfileobj(inp,dest)
        h=hashlib.sha256();n=0
        with gzip.open(packed,'rb') as inp:
            while b:=inp.read(1048576):h.update(b);n+=len(b)
        if h.hexdigest()!=raw_hash or n!=raw_bytes:raise ValueError('native final gzip roundtrip failed')
        src.unlink();retention.append(dict(raw_path=name,raw_sha256=raw_hash,raw_bytes=raw_bytes,gzip_path=str(packed.relative_to(args.output)),gzip_sha256=digest(packed),gzip_bytes=packed.stat().st_size))
    verify_observations(args.observed_receipt,args.geometry)
    complete_manifest(args.output,dict(stage='acquisition',geometry=args.geometry,block=args.block,enzyme=enzyme,run_id=run_id,physical_chains=selected.height,physical_H=selected.filter(pl.col('locus')=='IGH').height,physical_L=selected.filter(pl.col('locus')!='IGH').height,heavy_intensity=10000,absolute_heavy_band_recovery='unchanged heavy amount assumption; not measured gel recovery',candidate_count=inference['candidates'].height,background_records=inference['candidates'].height-12,query_records=12,retention=retention,retirement=retirement,command=command),inputs)


def _class_ids(sequences):
    from vdjaqp.match import ms_collapse
    return ['mass_'+hashlib.sha256(s.encode()).hexdigest() for s in ms_collapse(sequences)]



def experimental_readouts(panel,physical,queries,hits,blocks,native_by_block,targets):
    """Reporting only: reuse full-protein hits with supplied-query target masks."""
    import polars as pl
    from vdjaqp.airr_ms_metrics import supplied_reference_regions,variable_sequences,variable_only_placements,species_metrics
    from vdjaqp.abundance_inference import REGIONS,target_vector,distance
    observed=target_vector(targets)
    names=physical.filter(pl.col('locus')=='IGH').select('chain_id','query_clone_id',pl.col('reference').str.replace(r'^ref0*','ref').alias('reference'))
    if names.height!=12 or set(names['reference'])!={f'ref{i}' for i in range(1,13)}:raise ValueError('experimental reference mapping differs')
    supplied=panel.join(names.select(pl.col('query_clone_id').alias('clone_id'),'reference'),on='clone_id',validate='1:1').drop('clone_id').rename({'reference':'clone_id'})
    if supplied.height!=12:raise ValueError('supplied query panel population differs')
    regions=supplied_reference_regions(supplied);variables=variable_sequences(supplied)
    raw_to_species=queries.select('peptide_id',pl.col('seq').str.replace_all('I','L',literal=True).alias('seq'))
    species_queries=raw_to_species.select('seq').unique().sort('seq').with_row_index('_q').with_columns(pl.format('species_{}',pl.col('_q')).alias('_species_id')).drop('_q')
    remap=raw_to_species.join(species_queries,on='seq',validate='m:1')
    placements=(hits.join(names.select('chain_id','reference'),on='chain_id',how='inner').drop('chain_id').rename({'reference':'chain_id'})
        .join(remap.select('peptide_id','_species_id'),on='peptide_id').drop('peptide_id').rename({'_species_id':'peptide_id'}).select('chain_id','peptide_id','start','end').unique())
    placements=variable_only_placements(placements,variables.rename({'clone_id':'chain_id'}))
    species_queries=species_queries.rename({'_species_id':'peptide_id'})
    frames={};coverages={};inventories={};sim=[]
    views=[(b,blocks[b],len(native_by_block[b]))for b in BLOCKS]
    views.extend((condition,pl.concat([blocks[tryp],blocks[BLOCKS[2]]],how='vertical_relaxed'),len(native_by_block[tryp])+len(native_by_block[BLOCKS[2]])) for condition,tryp in [('solution_pooled',BLOCKS[0]),('heavy_band_pooled',BLOCKS[1])])
    displays={BLOCKS[0]:[('solution_mixture','trypsin')],BLOCKS[1]:[('gel_mixture','trypsin')],BLOCKS[2]:[('solution_mixture','chymotrypsin'),('gel_mixture','chymotrypsin')], 'solution_pooled':[('solution_mixture','pooled')],'heavy_band_pooled':[('gel_mixture','pooled')]}
    for label,calls,native_groups in views:
        inventory,coverage,pooled=species_metrics(calls,species_queries,placements,regions)
        inventory=pl.concat([inventory,pl.DataFrame(dict(quantity=['native_ms2_groups'],count=[native_groups],units=['spectrum_groups'],n_observations=[calls.height]))],how='vertical_relaxed')
        frames[label+'_query_species']=inventory;frames[label+'_supplied_query_coverage']=coverage;frames[label+'_supplied_query_coverage_summary']=pooled
        coverages[label]=coverage;inventories[label]=inventory
        for condition,enzyme in displays[label]:
            for row in pl.concat([coverage.select(pl.col('chain_id').alias('reference'),'region','covered_residues','region_length','coverage'),pooled.select(pl.lit('all').alias('reference'),'region','covered_residues','region_length','coverage')],how='vertical_relaxed').iter_rows(named=True):
                sim.append(dict(condition=condition,enzyme=enzyme,reference=row['reference'],region=row['region'],quantity='residue_union_coverage',simulated_value=row['coverage'],simulated_numerator=row['covered_residues'],simulated_denominator=row['region_length'],readout=label))
            value=inventory.filter(pl.col('quantity')=='compatible_species')['count'].item()
            sim.append(dict(condition=condition,enzyme=enzyme,reference='all',region=None,quantity='compatible_called_species',simulated_value=float(value),simulated_numerator=value,simulated_denominator=None,readout=label))
    summary=dict(yield_counts=[],coverage=[],denominators=[])
    for label in BLOCKS:
        summary['yield_counts'].append(inventories[label].filter(pl.col('quantity')=='compatible_species')['count'].item())
        lookup={(r['chain_id'],r['region']):r for r in coverages[label].to_dicts()}
        summary['coverage'].append([[lookup[(f'ref{i}',r)]['coverage']for r in REGIONS]for i in range(1,13)])
        summary['denominators'].append([[lookup[(f'ref{i}',r)]['region_length']for r in REGIONS]for i in range(1,13)])
    keys=['condition','enzyme','reference','region','quantity']
    real=targets.filter(pl.col('condition').is_in(['solution_mixture','gel_mixture']) & pl.col('quantity').is_in(['residue_union_coverage','compatible_called_species'])).select(*keys,'units',pl.col('value').alias('real_value'),pl.col('numerator').alias('real_numerator'),pl.col('denominator').alias('real_denominator'))
    predicted=pl.DataFrame(sim)
    if real.unique(keys).height!=real.height or predicted.unique(keys).height!=predicted.height:raise ValueError('duplicate experimental comparison key')
    residual=real.join(predicted,on=keys,how='full',coalesce=True,nulls_equal=True)
    if residual.filter(pl.col('real_value').is_null()|pl.col('simulated_value').is_null()|((pl.col('quantity')=='residue_union_coverage') & (pl.col('real_denominator')!=pl.col('simulated_denominator')))).height:raise ValueError('experimental residual keys/denominators differ')
    frames['experimental_target_residuals']=residual.with_columns((pl.col('simulated_value')-pl.col('real_value')).alias('signed_residual'),(pl.col('simulated_value')-pl.col('real_value')).abs().alias('absolute_residual')).sort(keys)
    return frames,dict(distance=distance(summary,observed),coordinates=255,summary=summary,observed=observed,unique_acquisitions=3,pooled_conditions=2,shared_chymotrypsin_selection_count=1,recovery_used_for_selection=False,scope='exact supplied-query heavy-variable compatibility, not native physical AA-inferred regions or total identified-species accuracy',attrition='noncanonical_rows and short_rows are overlapping diagnostic flags, not summed sequential losses; no PEAKS Area/ppm/ALC filter applied to native search scores')


def raw_background(args):
    import polars as pl
    v=json.loads(args.background_validation_receipt.read_text());raw,n=raw_truth_digest(args.background_truth_gzip)
    if v['status']!='complete' or (raw,n)!=(v['truth']['raw_sha256'],v['truth']['raw_bytes']) or digest(args.background_truth_gzip)!=v['truth']['gzip_sha256']:raise ValueError('raw background truth identity mismatch')
    with gzip.open(args.background_truth_gzip,'rt') as f:truth=json.load(f)
    if len(truth['chains'])!=v['observed_chains'] or truth['heavy_recombination_model']['sha256']!=v['model_sha256']:raise ValueError('raw background population/model mismatch')
    rows=[dict(chain_id=c['chain_id'],locus='IGH',variable_aa=c['variable_aa'],variable_length=len(c['variable_aa']),protein_aa=c['protein_aa']) for c in truth['chains'] if c['locus']=='IGH']
    if len(rows)!=v['locus_counts']['IGH']:raise ValueError('raw heavy population mismatch')
    return pl.DataFrame(rows)


def evaluate(args):
    import polars as pl
    import numpy as np
    from vdjaqp.readout import read_idxml
    from vdjaqp.simulation_metrics import exact_peptide_placements,identity_metrics
    from vdjaqp.pipeline import Evidence
    from vdjaqp.errormodel import ErrorModel,mass_prior,_E_DB_DEFAULT
    from vdjaqp.data import to_gene_level
    from vdjaqp.bench import classification_metrics,rank_targets
    from vdjaqp.airr_capture import capture_tables
    from scripts.smoke_metrics import aggregate_coverage,aggregate_gene_coverage
    verify_manifest(args.reference_regions)
    if json.loads(args.reference_regions.read_text())['packet_manifest_sha256']!=digest(args.packet_manifest):raise ValueError('region/packet identity mismatch')
    _,_,observed=verify_observations(args.observed_receipt,args.geometry)
    inference,_=inference_inputs(observed,args.packet_manifest);candidates=to_gene_level(inference['candidates'])
    physical=pl.read_parquet(args.reference_regions.parent/'physical_metadata.parquet');regions=pl.read_parquet(args.reference_regions.parent/'physical_regions.parquet')
    background=raw_background(args)
    all_physical=pl.concat([background,physical.select(background.columns)],how='vertical_relaxed')
    if all_physical['chain_id'].n_unique()!=all_physical.height:raise ValueError('evaluation chain identity collision')
    all_physical=all_physical.with_columns(pl.Series('variable_mass_class_id',_class_ids(all_physical['variable_aa'])))
    physical_ids=set(physical['chain_id']);heavy_ids=set(physical.filter(pl.col('locus')=='IGH')['chain_id'])
    target_records=set(inference['reference_panel']['sequence_id'])
    # AIRR-compatible reference_panel keeps source sequence_id, candidate rows keep clone_id.
    if not target_records<=set(candidates['clone_id']):raise ValueError('canonical reference identities differ')
    blocks={};native_by_block={};receipts={};inputs=input_pins([args.packet_manifest,args.reference_regions,args.observed_receipt,args.background_truth_gzip,args.background_validation_receipt,*args.acquisition_manifests,*([args.experimental_targets] if getattr(args,'experimental_targets',None) else [])])
    for path in args.acquisition_manifests:
        m=verify_manifest(path)
        if m['stage']!='acquisition' or m['geometry']!=args.geometry or m['block'] in blocks:raise ValueError('wrong/duplicate acquisition identity')
        _,heavy_only=acquisition_scope(m['block'])
        if m['physical_H']!=physical.filter(pl.col('locus')=='IGH').height or m['physical_L']!=(0 if heavy_only else physical.filter(pl.col('locus')!='IGH').height):raise ValueError('acquisition physical preparation mismatch')
        for p in [args.packet_manifest,args.reference_regions,args.observed_receipt]:
            if m['input_sha256'].get(str(p.resolve()))!=digest(p):raise ValueError('acquisition source mismatch')
        with gzip.open(path.parent/'native/psms.idXML.gz','rb') as f:blocks[m['block']]=read_idxml(f,run_id=m['run_id'])
        with gzip.open(path.parent/'native/ms2_truth.json.gz','rt') as f:native_by_block[m['block']]=json.load(f)
        receipts[m['block']]=m
    if set(blocks)!=set(BLOCKS):raise ValueError('exactly the three declared unique acquisitions are required')
    # ONE exact evaluation batch, with all distinct candidate/call/native peptide AA.
    candidate_aliases=candidates.select('clone_id','seq')
    sequences=pl.concat([candidate_aliases.select('seq'),*[f.select('seq') for f in blocks.values()],pl.DataFrame({'seq':[c['peptide_aa'] for events in native_by_block.values() for e in events for c in e['contributors']]},schema={'seq':pl.String})]).unique().sort('seq')
    valid=pl.col('seq').str.contains('^[ACDEFGHIKLMNPQRSTVWY]{6,}$')
    excluded=sequences.filter(~valid.fill_null(False));queries=sequences.filter(valid).with_row_index('_q').with_columns(pl.format('eval_{}',pl.col('_q')).alias('peptide_id')).drop('_q')
    started=time.monotonic();hits=exact_peptide_placements(all_physical.select('chain_id','protein_aa'),queries,threads=args.cpus);search_seconds=time.monotonic()-started
    if hits.height>args.max_placements:raise ValueError('native exact placement count exceeds explicit bound')
    aliases=candidate_aliases.join(queries,on='seq',how='inner')
    heavy=all_physical.filter(pl.col('locus')=='IGH')
    full=(hits.join(heavy.select('chain_id','variable_mass_class_id','variable_length'),on='chain_id')
          .filter((pl.col('start')==0)&(pl.col('end')>=pl.col('variable_length'))).select('peptide_id','variable_mass_class_id').unique())
    estimate=full.join(aliases.group_by('peptide_id').len(),on='peptide_id').select(pl.col('len').sum()).item() or 0
    if estimate>args.max_placements:raise ValueError('full-class alias expansion exceeds explicit bound')
    full_pairs=full.join(aliases.select('clone_id','peptide_id'),on='peptide_id').select('clone_id','variable_mass_class_id').unique()
    universe=heavy.select('variable_mass_class_id',pl.col('chain_id').is_in(heavy_ids).alias('physical_positive'))
    reference_hits=hits.filter(pl.col('chain_id').is_in(physical_ids))
    association_size=reference_hits.filter(pl.col('chain_id').is_in(heavy_ids)).join(aliases.group_by('peptide_id').len(),on='peptide_id').select(pl.col('len').sum()).item() or 0
    if association_size>args.max_placements:raise ValueError('physical candidate association expansion exceeds explicit bound')
    associations=reference_hits.filter(pl.col('chain_id').is_in(heavy_ids)).join(aliases.select('clone_id','peptide_id'),on='peptide_id').select('clone_id','chain_id','start','end').unique()
    args.output.mkdir(parents=True)
    provenance_columns=['clone_id','background_member','sequence_origin','coverage_kind','source_observation_id','observed_read_id','parent_clone_id','native_clone_row','duplicate_count','native_junction_nt','native_junction_aa']
    write_frames(args.output,{'candidate_source_provenance':candidates.select([c for c in provenance_columns if c in candidates.columns]),'evaluation_queries':queries,'excluded_short_queries':excluded,'exact_eval_placements':hits,'full_class_associations':full_pairs,'physical_reference_associations':associations})
    if getattr(args,'experimental_targets',None):
        _,panel,_=packet_frames(args.packet_manifest)
        comparison_frames,comparison=experimental_readouts(panel,physical,queries,hits,blocks,native_by_block,pl.read_csv(args.experimental_targets,separator='\t'))
        write_frames(args.output,comparison_frames)
        comparison_frames['experimental_target_residuals'].write_csv(args.output/'experimental_target_residuals.tsv',separator='\t',line_terminator='\n')
        comparison['experimental_targets_sha256']=digest(args.experimental_targets)
        (args.output/'experimental_comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
    results={};inventory=[]
    for block in BLOCKS:
        for row in physical.select('chain_id','locus').iter_rows(named=True):inventory.append(dict(block=block,**row,assayed=(block==BLOCKS[0] or row['locus']=='IGH')))
    write_frames(args.output,{'physical_assay_inventory':pl.DataFrame(inventory)})
    for condition,tryp in [('solution',BLOCKS[0]),('heavy_band',BLOCKS[1])]:
        calls=pl.concat([blocks[tryp],blocks[BLOCKS[2]]],how='vertical_relaxed')
        model=ErrorModel(np.array([0.,1.]),np.array([.1,.01]),mass_prior(),_E_DB_DEFAULT)
        # Only observed+supplied candidate rows and raw calls enter inference.
        evidence=Evidence.from_frames(candidates,calls,error_model=model);scores=evidence.score()
        if scores.height!=candidates.height or not scores['score'].is_finite().all():raise ValueError('invalid complete score population')
        ranks=rank_targets(scores,'score',set(candidates['clone_id']))
        identity,called,recovery=identity_metrics(physical.filter(pl.col('locus')=='IGH'),candidates,ranks,associations)
        classes,mapping,metrics=class_benchmark(universe,full_pairs,scores)
        metrics['canonical_query_row_retrieval']=classification_metrics(scores,target_records)
        metrics['scope']='canonical query ROW retrieval versus raw-generated+physical-reference full-variable I/L CLASS presence; recovery at20, no calibrated posterior claim'
        metrics['shared_chymotrypsin']='one actual acquisition contributes once to each condition union; unions are not independent replicates'
        results[condition]=metrics
        frames={'scores':scores,'ranks':ranks,'physical_identity_population':identity,'physical_candidate_identity':called,'physical_recovery':recovery,'class_scores':classes,'candidate_class_mapping':mapping}
        for name,frame in capture_tables(evidence,scores,target_records).items():frames['capture_'+name]=frame
        # Physical AA-derived coverage with explicit unassayed L status in heavy band.
        levels=['native_ms2_contributor_compatibility','identified_peptide_compatibility','accepted_heavy_candidate_exact_compatibility']
        selected=physical if condition=='solution' else physical.filter(pl.col('locus')=='IGH');selected_ids=set(selected['chain_id'])
        ref_hits=reference_hits.filter(pl.col('chain_id').is_in(selected_ids))
        identified_q=queries.join(calls.select('seq').unique(),on='seq',how='semi').select('peptide_id')
        called_q=aliases.filter(pl.col('clone_id').is_in(scores.filter(pl.col('confirmed_nats')>=20)['clone_id'].implode())).select('peptide_id').unique()
        links=[dict(seq=c['peptide_aa'],chain_id=protein) for block in [tryp,BLOCKS[2]] for event in native_by_block[block] for c in event['contributors'] for protein in c['protein_accessions']]
        native_links=pl.DataFrame(links,schema={'seq':pl.String,'chain_id':pl.String}).unique().join(queries,on='seq').select('peptide_id','chain_id')
        placed=pl.concat([ref_hits.join(native_links,on=['peptide_id','chain_id'],how='semi').select('chain_id','start','end').unique().with_columns(pl.lit(levels[0]).alias('evidence_level')),
                         ref_hits.join(identified_q,on='peptide_id',how='semi').select('chain_id','start','end').unique().with_columns(pl.lit(levels[1]).alias('evidence_level')),
                         ref_hits.join(called_q,on='peptide_id',how='semi').select('chain_id','start','end').unique().with_columns(pl.lit(levels[2]).alias('evidence_level'))])
        coverage,cost=bounded_coverage(regions,placed,levels,max_slots=args.max_expanded_slots)
        if condition=='heavy_band':
            absent=pl.col('locus')!='IGH'
            coverage=coverage.with_columns(*[pl.when(absent).then(None).otherwise(pl.col(k)).alias(k) for k in ['region_length','covered_residues','coverage']],*[pl.when(absent).then(pl.lit('not_assayed')).otherwise(pl.col(k)).alias(k) for k in ['status','region_status']])
        frames.update(physical_region_coverage=coverage,physical_gene_coverage=aggregate_gene_coverage(coverage,physical),physical_locus_coverage=aggregate_coverage(coverage,['locus','region','evidence_level']),coverage_expansion=cost)
        write_frames(args.output,{condition+'_'+n:f for n,f in frames.items()})
        del evidence
    (args.output/'metrics.json').write_text(json.dumps(results,indent=2)+'\n')
    verify_observations(args.observed_receipt,args.geometry)
    complete_manifest(args.output,dict(stage='pooled_evaluation',geometry=args.geometry,unique_native_acquisitions=3,pooled_conditions=2,shared_chymo_count_once=True,actual_candidates=candidates.height,raw_background_heavy_records=background.height,evaluation_query_count=queries.height,exact_eval_batches=1,exact_placement_rows=hits.height,native_exact_eval_seconds=search_seconds,full_class_alias_estimate=estimate,physical_association_estimate=association_size,physical_unassayed_L='heavy-band region rows null/not_assayed; no zero-coverage claim'),inputs)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='stage',required=True)
    for stage in ['qualify-regions','acquire','evaluate']:
        p=sub.add_parser(stage);p.add_argument('--packet-manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--cpus',type=int,default=16)
        if stage!='qualify-regions':
            p.add_argument('--observed-receipt',type=Path,required=True);p.add_argument('--reference-regions',type=Path,required=True);p.add_argument('--geometry',choices=['amplicon','rnaseq'],required=True)
        if stage=='acquire':
            p.add_argument('--block',choices=BLOCKS,required=True);p.add_argument('--mssimulator',type=Path,required=True);p.add_argument('--stage-timeout',type=float,default=900.)
        if stage=='evaluate':
            p.add_argument('--experimental-targets',type=Path,help='optional recorded experimental yield/coverage TSV; reporting only')
            p.add_argument('--acquisition-manifests',type=Path,nargs=3,required=True);p.add_argument('--background-truth-gzip',type=Path,required=True);p.add_argument('--background-validation-receipt',type=Path,required=True)
            p.add_argument('--max-placements',type=int,default=5000000);p.add_argument('--max-expanded-slots',type=int,default=2000000)
    args=parser.parse_args(argv)
    if args.output.exists():raise FileExistsError('fresh output required')
    if getattr(args,'stage_timeout',1)<=0 or getattr(args,'max_placements',1)<1 or getattr(args,'max_expanded_slots',1)<1:raise ValueError('positive stage and memory bounds required')
    from vdjaqp import resolve_cpu_budget
    args.cpus=resolve_cpu_budget(args.cpus)
    import os
    os.environ['ARDA_NO_AUTO_FETCH']='1'
    return {'qualify-regions':qualify_regions,'acquire':acquire,'evaluate':evaluate}[args.stage](args)


if __name__=='__main__':main()
