"""Describe one validated generated substrate across heavy, kappa and lambda.

Raw observed chains only; no generation, annotation, fitting, inferred ancestry,
mutation reconstruction, assay calculation or experimental abundance inference.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))


def digest(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def summarize(truth, *, expected_chains=66080):
    """Project native fields once, then use vectorized table summaries."""
    import polars as pl
    from vdjaqp.data import to_gene_level
    from scripts.summarize_simulation_study import gene_frequencies
    chains=truth['chains'];ancestors=truth.get('latent_ancestors',[])
    if len(chains)!=expected_chains:raise ValueError('observed chain count differs from declared cohort')
    rows=[]
    def optional_count(value):
        if value is not None and (type(value) is not int or value<0):
            raise ValueError('reported count must be a nonnegative integer or missing')
        return value
    for chain in chains:
        a=chain['airr_truth']
        if chain.get('observed',True) is not True or a.get('productive') is not True:
            raise ValueError('observed chains must be productive; latent ancestors are excluded')
        required=['chain_id','antibody_id','locus','variable_nt','coding_variable_nt','variable_aa','protein_aa']
        if any(not isinstance(chain.get(k),str) or not chain[k] for k in required):
            raise ValueError('missing native chain identity/sequence fields')
        jnt,jaa=a.get('junction'),a.get('junction_aa')
        if (chain['locus'] not in ['IGH','IGK','IGL'] or not isinstance(jnt,str) or not isinstance(jaa,str)
                or not jaa.startswith('C') or jaa[-1] not in 'FW' or len(jnt)!=3*len(jaa)):
            raise ValueError('native locus or anchor-inclusive junction schema is invalid')
        row={k:chain.get(k) for k in ['chain_id','antibody_id','locus','family_id','component','constant_nt_gene']}
        row.update(v_call=a.get('v_call'),j_call=a.get('j_call'),junction_nt=jnt,junction_aa=jaa,
            variable_aa=chain['variable_aa'],variable_nt_length=len(chain['variable_nt']),
            coding_variable_nt_length=len(chain['coding_variable_nt']),variable_aa_length=len(chain['variable_aa']),
            protein_aa_length=len(chain['protein_aa']),anchored_junction_nt_length=len(jnt),
            anchored_junction_aa_length=len(jaa),native_n_mutations=optional_count(a.get('n_mutations')))
        for relation in ['root','branch']:
            distance=chain.get(relation+'_distance')
            if distance is not None and not isinstance(distance,dict):raise ValueError('distance descriptor must be native fields or missing')
            distance=distance or {}
            for name in ['substitutions','denominator_nt']:
                row[relation+'_'+name]=optional_count(distance.get(name))
            known=[row[relation+'_'+name] for name in ['substitutions','denominator_nt']]
            if all(x is not None for x in known) and known[0]>known[1]:raise ValueError('endpoint substitutions exceed declared NT denominator')
        rows.append(row)
    strings=['chain_id','antibody_id','locus','family_id','component','constant_nt_gene','v_call','j_call','junction_nt','junction_aa','variable_aa']
    frame=pl.DataFrame(rows,schema_overrides={**{k:pl.String for k in strings},**{k:pl.Int64 for k in rows[0] if k not in strings}})
    if frame['chain_id'].n_unique()!=frame.height:raise ValueError('duplicate observed chain ID')
    observed_ids=set(frame['chain_id'])
    if any(a.get('observed') is not False or a.get('chain_id') in observed_ids for a in ancestors):
        raise ValueError('latent ancestor membership conflicts with observed population')
    pairs=frame.group_by('antibody_id').agg(pl.len().alias('chains'),(pl.col('locus')=='IGH').sum().alias('heavy'))
    if pairs.filter((pl.col('chains')!=2)|(pl.col('heavy')!=1)).height:raise ValueError('each antibody requires one heavy and one light chain')
    counts=dict(frame.group_by('locus').len().iter_rows())
    if expected_chains==66080 and counts!={'IGH':33040,'IGK':16520,'IGL':16520}:raise ValueError('matched cohort locus counts differ')
    frame=frame.with_columns(pl.col('component').fill_null('missing'),pl.col('junction_aa').str.replace_all('I','L',literal=True).alias('junction_IL'),pl.col('variable_aa').str.replace_all('I','L',literal=True).alias('variable_IL'))
    population=frame.group_by('locus','component').agg(pl.len().alias('observed_chain_records'),pl.col('antibody_id').n_unique().alias('antibody_records'),pl.col('family_id').null_count().alias('missing_family_records'))
    calls=frame.select('locus',pl.col('v_call').alias('v_name'),pl.col('j_call').alias('j_name'),pl.col('constant_nt_gene').alias('c_name'))
    fields={'V':'v_name','J':'j_name','C':'c_name'}
    genes=pl.concat([gene_frequencies(calls,fields,'observed_generated_chains','native_GenAIRR_allele_calls'),
        gene_frequencies(to_gene_level(calls),fields,'observed_generated_chains','native_GenAIRR_gene_calls')])
    lengths=[]
    for field in ['variable_nt_length','coding_variable_nt_length','variable_aa_length','protein_aa_length','anchored_junction_nt_length','anchored_junction_aa_length']:
        counts=frame.group_by('locus',pl.col(field).alias('value')).len(name='count')
        lengths.append(counts.with_columns(pl.lit(field).alias('quantity'),pl.lit('nt' if '_nt_' in field else 'aa').alias('unit')))
    lengths=pl.concat(lengths).join(frame.group_by('locus').len(name='denominator'),on='locus').with_columns((pl.col('count')/pl.col('denominator')).alias('fraction'))
    identities=frame.group_by('locus').agg(pl.len().alias('chain_records'),
        pl.col('junction_nt').n_unique().alias('anchored_junction_NT_classes'),
        pl.col('junction_aa').n_unique().alias('anchored_junction_AA_classes'),
        pl.col('junction_IL').n_unique().alias('anchored_junction_AA_IL_classes'),
        pl.col('variable_IL').n_unique().alias('variable_AA_IL_classes'))
    families=frame.filter(pl.col('family_id').is_not_null()).group_by('locus','component','family_id').len(name='observed_chains_per_family')
    family_sizes=families.group_by('locus','component','observed_chains_per_family').len(name='family_count').join(
        families.group_by('locus','component').len(name='known_family_denominator'),on=['locus','component'])
    shm=[];hist=[]
    for relation in ['root','branch']:
        substitution,denominator=relation+'_substitutions',relation+'_denominator_nt'
        valid=pl.col(substitution).is_not_null()&pl.col(denominator).is_not_null()
        part=frame.group_by('locus','component').agg(pl.len().alias('chain_records'),valid.sum().alias('available_endpoint_records'),
            pl.col(substitution).filter(valid).sum().alias('substitutions'),pl.col(denominator).filter(valid).sum().alias('denominator_nt'))
        shm.append(part.with_columns(pl.lit(relation).alias('endpoint'),
            (pl.col('chain_records')-pl.col('available_endpoint_records')).alias('missing_endpoint_records'),
            *[pl.when(pl.col('available_endpoint_records')>0).then(pl.col(k)).otherwise(None).alias(k) for k in ['substitutions','denominator_nt']])
            .with_columns(pl.when(pl.col('denominator_nt')>0).then(pl.col('substitutions')/pl.col('denominator_nt')).otherwise(None).alias('endpoint_substitution_fraction')))
    for field in ['root_substitutions','branch_substitutions','native_n_mutations']:
        part=frame.group_by('locus','component',pl.col(field).alias('value')).len(name='count')
        hist.append(part.with_columns(pl.lit(field).alias('quantity'),pl.col('value').is_null().alias('missing_value')))
    hist=pl.concat(hist).join(frame.group_by('locus','component').len(name='denominator'),on=['locus','component'])
    return dict(population_counts=population,gene_usage=genes,length_distributions=lengths,identity_counts=identities,
        family_sizes=family_sizes,shm_descriptors=pl.concat(shm),shm_distributions=hist)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    for flag in ['truth-gzip','validation-receipt','output']:parser.add_argument('--'+flag,type=Path,required=True)
    parser.add_argument('--cpus',type=int,default=16)
    args=parser.parse_args(argv)
    if args.output.exists() or args.output.is_symlink():raise FileExistsError('fresh output required')
    from vdjaqp import resolve_cpu_budget
    cpus=resolve_cpu_budget(args.cpus)
    from scripts.export_mock_observation_metrics import raw_truth_digest
    validation=json.loads(args.validation_receipt.read_text())
    if validation.get('status')!='complete' or validation['observed_chains']!=66080:raise ValueError('completed full-cohort validation receipt required')
    pins={str(p.resolve()):digest(p) for p in [args.truth_gzip,args.validation_receipt]}
    raw,rawbytes=raw_truth_digest(args.truth_gzip)
    if (raw,rawbytes)!=(validation['truth']['raw_sha256'],validation['truth']['raw_bytes']) or pins[str(args.truth_gzip.resolve())]!=validation['truth']['gzip_sha256']:
        raise ValueError('validated native truth raw/gzip identity differs')
    sources=[Path(__file__),Path(__file__).with_name('summarize_simulation_study.py'),Path(__file__).with_name('export_mock_observation_metrics.py'),Path(__file__).parents[1]/'vdjaqp/data.py']
    source_pins={str(p.resolve()):digest(p) for p in sources}
    with gzip.open(args.truth_gzip,'rt') as stream:truth=json.load(stream)
    frames=summarize(truth)
    args.output.mkdir(parents=True)
    for name,frame in frames.items():frame.sort(frame.columns).write_csv(args.output/(name+'.tsv'),separator='\t',line_terminator='\n')
    if any(digest(Path(p))!=h for p,h in {**pins,**source_pins}.items()):raise ValueError('source/input drift')
    manifest=dict(status='complete',cpu_budget=cpus,observed_chains=66080,excluded_latent_chains=len(truth.get('latent_ancestors',[])),
        input_sha256=pins,source_sha256=source_pins,artifact_sha256={p.name:digest(p) for p in args.output.iterdir()},
        table_rows={name:frame.height for name,frame in frames.items()},raw_truth_sha256=raw,
        source_model_sha256=truth.get('heavy_recombination_model',{}).get('sha256'),
        scope='Descriptive native genotype calls; all-locus sequence classes are not unique B cells. Light composition and constants are model assumptions, not empirical fits. Existing heavy calibration/validation metrics remain authoritative.',
        family_scope='Declared family IDs/components and reported parent/root endpoint counts only; missing fields remain unavailable, no ancestry or mutation inference.',
        abundance_scope=truth.get('abundance_scope'),lambda_scope=truth.get('lambda_scope'))
    (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


if __name__=='__main__':main()
