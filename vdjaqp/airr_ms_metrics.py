"""Evaluation-only supplied AIRR boundaries and exact called-species readouts."""
import polars as pl
from .simulation_metrics import region_coverage


def variable_sequences(candidates):
    """Remove only the supplied contiguous constant-labelled suffix, never fill bases."""
    rows=[]
    for row in candidates.select('clone_id','seq','region_ids').iter_rows(named=True):
        seq,labels=row['seq'],row['region_ids']
        if len(seq)!=len(labels) or any(type(z) is not int or z not in range(9) for z in labels):
            raise ValueError('invalid supplied residue labels')
        n=labels.count(8)
        end=len(seq)-n
        if (n and labels[end:]!=[8]*n) or 8 in labels[:end] or end==0:
            raise ValueError('constant labels must be a contiguous suffix of a nonempty variable')
        rows.append((row['clone_id'],seq[:end],end))
    return pl.DataFrame(rows,schema={'clone_id':pl.String,'variable_aa':pl.String,'variable_length':pl.Int64},orient='row')


def supplied_reference_regions(references):
    """Reconcile supplied FR/CDR boundaries with anchor-excluded CDR3, as in real targets."""
    variables=variable_sequences(references)
    lengths=dict(zip(variables['clone_id'],variables['variable_length']))
    rows=[]
    for x in references.iter_rows(named=True):
        seq=x['seq'][:lengths[x['clone_id']]]; labels=x['region_ids']
        b,e=x['cdr3_start'],x['cdr3_end'];junction=x['cdr3_seq']
        if not junction or seq[b:e]!=junction or seq.count(junction)!=1 or junction[0]!='C' or junction[-1] not in 'FW':
            raise ValueError('invalid supplied anchor-inclusive junction')
        if [i for i,z in enumerate(labels) if z in (5,6)]!=list(range(b,e)):
            raise ValueError('supplied junction labels disagree')
        regions=[]
        for name,label in [('FR1',0),('CDR1',1),('FR2',2),('CDR2',3),('FR3',4),('FR4',7)]:
            positions=[i for i,z in enumerate(labels) if z==label]
            if not positions:
                regions.append((name,None,None));continue
            lo,hi=positions[0],positions[-1]+1
            if positions!=list(range(lo,hi)):raise ValueError('noncontiguous supplied region')
            if label==4:
                if hi!=b:raise ValueError('FR3 junction boundary conflict')
                hi=b+1
            if label==7:
                if lo!=e:raise ValueError('FR4 junction boundary conflict')
                lo=e-1
            regions.append((name,lo,hi))
        regions.extend([('CDR3',b+1,e-1),('junction',b,e),('variable',0,len(seq))])
        partition=sorted((lo,hi) for name,lo,hi in regions if name not in ('junction','variable') and lo is not None)
        if len(partition)==7 and (partition[0][0]!=0 or partition[-1][1]!=len(seq) or any(a[1]!=b[0] for a,b in zip(partition,partition[1:]))):
            raise ValueError('supplied partition does not cover variable')
        rows.extend((x['clone_id'],'IGH',name,lo,hi) for name,lo,hi in regions)
    return pl.DataFrame(rows,schema={'chain_id':pl.String,'locus':pl.String,'region':pl.String,'start':pl.Int64,'end':pl.Int64},orient='row')


def eligible_species(peptides):
    """Canonical >=6-AA I/L-equivalent distinct species, not PSM assignments."""
    return (peptides.filter(pl.col('seq').str.contains('^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)&
                            (pl.col('seq').str.len_chars()>=6))
            .select(pl.col('seq').str.replace_all('I','L',literal=True).alias('seq')).unique().sort('seq'))


def species_metrics(peptides, queries, placements, regions):
    """One supplied physical reference population; exact unions retain all denominators."""
    species=eligible_species(peptides).join(queries,on='seq',validate='1:1')
    hits=placements.join(species.select('peptide_id'),on='peptide_id',how='semi')
    coverage=region_coverage(regions,hits.with_columns(pl.lit('called_species').alias('evidence_level')),['called_species'])
    compatible=hits['peptide_id'].n_unique()
    canonical=peptides['seq'].str.contains('^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)
    counts={'psm_assignment_rows':peptides.height,'noncanonical_rows':int((~canonical).sum()),
            'short_rows':peptides.filter(pl.col('seq').str.len_chars()<6).height,
            'eligible_species':species.height,'compatible_species':compatible,
            'unmapped_species':species.height-compatible}
    inventory=pl.DataFrame({'quantity':list(counts),'count':list(counts.values()),
        'units':['rows' if key.endswith('rows') else 'species' for key in counts],
        'n_observations':[peptides.height]*len(counts)})
    pooled=coverage.group_by('region').agg(pl.len().alias('physical_reference_count'),
        pl.col('coverage').count().alias('defined_reference_count'),
        pl.col('covered_residues').sum().alias('covered_residues'),pl.col('region_length').sum().alias('region_length'))
    pooled=pooled.with_columns(*[pl.when(pl.col('defined_reference_count')>0).then(pl.col(key))
        .otherwise(None).alias(key) for key in ('covered_residues','region_length')])
    pooled=pooled.with_columns(pl.when(pl.col('defined_reference_count')>0)
        .then(pl.col('covered_residues')/pl.col('region_length')).otherwise(None).alias('coverage'))
    return inventory,coverage,pooled


def supplied_raw_regions(candidates):
    """Each supplied-label contiguous run, plus full/variable/stub views; no fill.

    Labels 5/6 retain producer meaning. They do not establish IMGT CDR3 or
    separate templated D from N additions. D sequence itself is templated.
    """
    variables=variable_sequences(candidates)
    lengths=dict(zip(variables['clone_id'],variables['variable_length']))
    names=['rawFR1','rawCDR1','rawFR2','rawCDR2','rawFR3','rawpart5','rawpart6','rawFR4','Cstub']
    rows=[]
    for x in candidates.iter_rows(named=True):
        key=x['clone_id'];end=lengths[key];labels=x['region_ids'];n=len(x['seq'])
        views=[('full_supplied','full_supplied',0,n,None),('variable','variable',0,end,None),
               ('Cstub','Cstub',end if end<n else None,n if end<n else None,x.get('c_name'))]
        for label,name in enumerate(names[:8]):
            positions=[i for i,z in enumerate(labels) if z==label]
            runs=[]
            for pos in positions:
                if runs and pos==runs[-1][1]:runs[-1][1]+=1
                else:runs.append([pos,pos+1])
            gene=x.get('v_name') if label<5 else x.get('j_name') if label==7 else None
            views.extend((f'{name}_run{i+1}',name,lo,hi,gene) for i,(lo,hi) in enumerate(runs))
            if not runs:views.append((name+'_run1',name,None,None,gene))
        rows.extend((key,'IGH',name,raw,lo,hi,gene.split('*')[0] if gene else None)
                    for name,raw,lo,hi,gene in views)
    return pl.DataFrame(rows,schema={'chain_id':pl.String,'locus':pl.String,'region':pl.String,
        'raw_region':pl.String,'start':pl.Int64,'end':pl.Int64,'gene_call':pl.String},orient='row')


def variable_only_placements(placements,metadata):
    """Exclude entire cross-boundary peptides; never clip them into variable hits."""
    return (placements.join(metadata.select('chain_id','variable_length'),on='chain_id',validate='m:1')
            .filter(pl.col('end')<=pl.col('variable_length')).drop('variable_length'))


def supplied_gene_coverage(coverage):
    """Pool residue unions per raw run and supplied gene, preserving missing genes."""
    return coverage.group_by(['locus','region','raw_region','gene_call','evidence_level']).agg(
        pl.col('chain_id').n_unique().alias('physical_chain_count'),
        pl.col('coverage').count().alias('defined_chain_count'),
        pl.col('covered_residues').sum().alias('covered_residues'),
        pl.col('region_length').sum().alias('residue_denominator')
    ).with_columns(*[pl.when(pl.col('defined_chain_count')>0).then(pl.col(key)).otherwise(None).alias(key)
                     for key in ('covered_residues','residue_denominator')]
    ).with_columns((pl.col('covered_residues')/pl.col('residue_denominator')).alias('pooled_residue_coverage'))
