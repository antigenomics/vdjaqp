"""Evaluation-only generated-chain metadata and annotation-derived residue regions.

No returned truth fields belong in inference inputs. Native ARDA NT coordinates
are 1-based inclusive; exported protein intervals are 0-based half-open.
"""
from __future__ import annotations

import polars as pl
from arda.refbuild.translate import translate

_PARTITION = (('fwr1', 'FR1'), ('cdr1', 'CDR1'), ('fwr2', 'FR2'),
              ('cdr2', 'CDR2'), ('fwr3', 'FR3'), ('cdr3', 'CDR3'), ('fwr4', 'FR4'))


def regions_from_annotations(chains, annotations):
    """Return ``(chain_metadata, regions)`` after validating one native NT batch.

    Missing annotation boundaries stay null. Generator junction placement defines
    overlapping views and validates, rather than repairs, inferred FR/CDR anchors.
    """
    chains = list(chains)
    annotations = list(annotations)
    ids = [c['chain_id'] for c in chains]
    ann_ids = [a['sequence_id'] for a in annotations]
    if len(set(ids)) != len(ids) or len(set(ann_ids)) != len(ann_ids) or set(ids) != set(ann_ids):
        raise ValueError('chain/annotation keys must be unique and equal')
    by_id = {a['sequence_id']: a for a in annotations}
    metadata, regions = [], []
    for c in chains:
        a = by_id[c['chain_id']]
        nt, var, protein = c['variable_nt'], c['variable_aa'], c['protein_aa']
        if c['locus'] not in ('IGH', 'IGK', 'IGL') or not var or set(protein) - set('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError('invalid chain locus/protein')
        if not protein.startswith(var) or len(protein) <= len(var):
            raise ValueError('physical protein must contain variable plus constant sequence')
        if translate(nt) != var or len(nt) not in range(3*len(var), 3*len(var)+3):
            raise ValueError('native coding variable translation/frame mismatch')
        if a.get('sequence') != nt or a.get('rev_comp') not in ('F', False, ''):
            raise ValueError('annotation must describe submitted forward coding variable')
        if c.get('rnaseq_nt'):
            rna = c['rnaseq_nt']
            if c.get('rna_model') == 'leader_imputed_v1':
                prefix, leader = c.get('rna_prefix_nt'), c.get('rna_leader_nt')
                start, initiation = c.get('rna_mature_nt_start'), c.get('rna_leader_nt_start')
                if (not isinstance(prefix, str) or not isinstance(leader, str) or not leader or
                        type(start) is not int or type(initiation) is not int or
                        start != len(prefix)+1 or initiation != len(prefix)-len(leader)+1 or
                        not prefix.endswith(leader) or not rna.startswith(prefix+nt) or
                        not leader.startswith('ATG') or len(leader)%3 or '*' in translate(leader) or
                        translate(rna[start-1:]) != protein+'*' or
                        translate(rna[initiation-1:]) != translate(leader)+protein+'*'):
                    raise ValueError('leader RNA/mature substrate coordinates or translation disagree')
            elif c.get('rna_model') not in (None, ''):
                raise ValueError('unknown RNA model')
            elif translate(rna) != protein+'*':
                raise ValueError('full native coding RNA does not translate to physical protein plus stop')
        genotype = c['airr_truth']
        junction = genotype['junction_aa'].upper()
        if not junction or var.count(junction) != 1 or junction[0] != 'C' or junction[-1] not in 'FW':
            raise ValueError('generator junction must have uniquely placed conserved anchors')
        jstart, jend = var.index(junction), var.index(junction)+len(junction)
        raw_junction = genotype.get('junction')
        raw_start, raw_end = genotype.get('junction_start'), genotype.get('junction_end')
        if (not isinstance(raw_junction, str) or not raw_junction
                or type(raw_start) is not int or type(raw_end) is not int
                or not 0 <= raw_start < raw_end <= len(nt)
                or raw_start % 3 or raw_end % 3
                or nt[raw_start:raw_end] != raw_junction.upper()
                or translate(raw_junction.upper()) != junction
                or (raw_start//3, raw_end//3) != (jstart, jend)):
            raise ValueError('native genotype junction coordinates/sequence/translation disagree')
        annotated_junction = a.get('junction_aa') or ''
        annotated_nt = a.get('junction') or ''
        concordance = ('missing' if not annotated_junction and not annotated_nt else
                       'conflict' if ((annotated_junction and annotated_junction != junction)
                                      or (annotated_nt and annotated_nt.upper() != raw_junction.upper())) else 'matches')
        meta = {key: c.get(key) for key in ('chain_id','locus','antibody_id','family_id','parent_id','constant_nt_gene','constant_accession','rna_pairs','secretion_factor','component','ancestor_id')}
        meta.update(variable_aa=var, protein_aa=protein, variable_length=len(var),
                    protein_length=len(protein), untranslated_terminal_nt=nt[3*len(var):],
                    v_call=genotype.get('v_call'), d_call=genotype.get('d_call'),
                    j_call=genotype.get('j_call'), n_mutations=genotype.get('n_mutations'),
                    n_observable_mutations=genotype.get('n_observable_mutations'),
                    generator_junction_nt=raw_junction, generator_junction_aa=junction,
                    generator_junction_start_nt=raw_start, generator_junction_end_nt=raw_end,
                    annotated_junction_nt=annotated_nt, annotated_junction_aa=annotated_junction,
                    junction_annotation_concordance=concordance,
                    annotated_locus=a.get('locus') or None, annotated_v_call=a.get('v_call') or None, annotated_j_call=a.get('j_call') or None)
        gene=c.get('canonical_v_gene'); status=c.get('canonical_v_match_status','not_evaluated')
        count=c.get('canonical_v_gene_count')
        valid = ((status=='unique_gene' and c['locus']=='IGH' and isinstance(gene,str) and gene.startswith('IGHV') and type(count) is int and count==1)
                 or (status=='ambiguous_genes' and c['locus']=='IGH' and gene is None and type(count) is int and count>1)
                 or (status=='no_exact_match' and c['locus']=='IGH' and gene is None and type(count) is int and count==0)
                 or (status=='not_evaluated' and gene is None and count is None)
                 or (status=='not_evaluated_locus' and c['locus']!='IGH' and gene is None and count is None))
        if not valid:raise ValueError('inconsistent canonical V mapping metadata')
        meta.update(canonical_v_gene=gene,canonical_v_match_status=status,canonical_v_gene_count=count)
        for label, field in [('branch','branch_distance'), ('ancestor','root_distance')]:
            distance = c.get(field) or {}
            meta.update({label+'_'+key: distance.get(source) for key,source in
                         [('substitutions','substitutions'),('denominator_nt','denominator_nt'),('fraction','fraction')]})
        metadata.append(meta)
        bounds, region_status = {}, {}
        for native, region in _PARTITION:
            start, end = a.get(native+'_start'), a.get(native+'_end')
            if start in (None, '') or end in (None, ''):
                bounds[region] = (None, None)
                region_status[region] = 'unavailable'
                continue
            start, end = int(start)-1, int(end)
            # Native FR4 can include the terminal split J-C codon. Only full
            # variable codons project to residues; retain its bases in metadata.
            if region == 'FR4' and end == len(nt):
                end = min(end, 3*len(var))
            if start < 0 or end <= start or end > len(var)*3:
                raise ValueError(f'invalid native coding region: {region}')
            if start % 3 or end % 3:
                bounds[region] = (None, None)
                region_status[region] = 'unavailable_non_codon_boundary'
                continue
            bounds[region] = (start//3, end//3)
            region_status[region] = 'available'
        known = [v for v in bounds.values() if v[0] is not None]
        if any(left[1] > right[0] for left, right in zip(known, known[1:])):
            raise ValueError('overlapping disjoint variable regions')
        for region, expected, endpoint in [('FR3',jstart+1,1), ('CDR3',jstart+1,0), ('CDR3',jend-1,1), ('FR4',jend-1,0)]:
            if bounds[region][endpoint] is not None and bounds[region][endpoint] != expected:
                bounds[region] = (None, None)
                region_status[region] = 'unavailable_junction_boundary_conflict'
        known = [v for v in bounds.values() if v[0] is not None]
        if len(known) == 7 and any(l[1] != r[0] for l,r in zip(known,known[1:])):
            raise ValueError('complete region partition must cover full variable')
        for region, endpoint, expected in [('FR1',0,0), ('FR4',1,len(var))]:
            if bounds[region][endpoint] is not None and bounds[region][endpoint] != expected:
                bounds[region] = (None, None)
                region_status[region] = 'unavailable_incomplete_variable_partition'
        known = [v for v in bounds.values() if v[0] is not None]
        if len(known) == 7 and (known[0][0] != 0 or known[-1][1] != len(var) or any(l[1] != r[0] for l,r in zip(known,known[1:]))):
            raise ValueError('complete region partition must cover full variable')
        for region, (start, end) in [*bounds.items(), ('JUNCTION',(jstart,jend)),
                                     ('GENERATOR_CDR3',(jstart+1,jend-1)), ('C',(len(var),len(protein)))]:
            regions.append(dict(chain_id=c['chain_id'],locus=c['locus'],region=region,start=start,end=end,
                                status=region_status.get(region, 'available'),
                                source='generator_junction' if region in ('JUNCTION','GENERATOR_CDR3') else 'physical_constant_boundary' if region=='C' else 'arda_nt_annotation',
                                gene_call=c.get('constant_nt_gene') if region=='C' else None))
        # GenAIRR exports native zero-based half-open segment coordinates,
        # unlike the one-based ARDA annotation endpoints above.
        for segment in ('V','D','J'):
            name=segment.lower()
            raw_start,raw_end=genotype.get(name+'_sequence_start'),genotype.get(name+'_sequence_end')
            start=end=None
            status='not_applicable' if segment=='D' and c['locus']!='IGH' else 'unavailable'
            if status!='not_applicable' and raw_start is not None and raw_end is not None:
                if not isinstance(raw_start,int) or not isinstance(raw_end,int) or not 0<=raw_start<=raw_end<=len(nt):
                    raise ValueError('invalid native genotype segment interval')
                lo,hi=(raw_start+2)//3,raw_end//3
                if hi>lo:
                    start,end,status=lo,hi,'available'
            regions.append(dict(chain_id=c['chain_id'],locus=c['locus'],region=segment,start=start,end=end,status=status,source='generator_fully_templated_codons',gene_call=genotype.get(name+'_call') or None))
    meta_schema = {'chain_id':pl.String,'locus':pl.String,'antibody_id':pl.String,'family_id':pl.String,'parent_id':pl.String,'constant_nt_gene':pl.String,'constant_accession':pl.String,'rna_pairs':pl.Int64,'secretion_factor':pl.Float64,'variable_aa':pl.String,'protein_aa':pl.String,'variable_length':pl.Int64,'protein_length':pl.Int64,'untranslated_terminal_nt':pl.String,'v_call':pl.String,'d_call':pl.String,'j_call':pl.String,'n_mutations':pl.Int64,'n_observable_mutations':pl.Int64,'annotated_locus':pl.String,'annotated_v_call':pl.String,'annotated_j_call':pl.String}
    meta_schema.update(canonical_v_gene=pl.String,canonical_v_match_status=pl.String,canonical_v_gene_count=pl.Int64)
    meta_schema.update(generator_junction_nt=pl.String,generator_junction_aa=pl.String,
                       generator_junction_start_nt=pl.Int64,generator_junction_end_nt=pl.Int64,
                       annotated_junction_nt=pl.String,annotated_junction_aa=pl.String,
                       junction_annotation_concordance=pl.String)
    meta_schema.update(component=pl.String,ancestor_id=pl.String,
                       branch_substitutions=pl.Int64,branch_denominator_nt=pl.Int64,branch_fraction=pl.Float64,
                       ancestor_substitutions=pl.Int64,ancestor_denominator_nt=pl.Int64,ancestor_fraction=pl.Float64)
    region_schema = {'chain_id':pl.String,'locus':pl.String,'region':pl.String,'start':pl.Int64,'end':pl.Int64,'status':pl.String,'source':pl.String,'gene_call':pl.String}
    return pl.DataFrame(metadata,schema=meta_schema), pl.DataFrame(regions,schema=region_schema)


def annotate_truth(chains, threads=None):
    """Annotate every complete generated variable in one evaluation-only NT batch."""
    from .resources import resolve_cpu_budget
    from arda.annotate.contig import reannotate_contigs
    chains = list(chains)
    budget = resolve_cpu_budget(threads)
    annotations = reannotate_contigs([(c['chain_id'],c['variable_nt']) for c in chains],
                                    seqtype='nt',threads=budget) if chains else []
    return regions_from_annotations(chains, annotations)
