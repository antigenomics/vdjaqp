"""Observation-only contig adapter; native weights and simulated truth are unused."""
from itertools import groupby
import polars as pl
from arda.refbuild.translate import translate, reverse_complement
from .match import MIN_QUERY_LEN

_CANDIDATE_SCHEMA={'clone_id':pl.String,'seq':pl.String,'nt_seq':pl.String,
 'cdr3_seq':pl.String,'cdr3_start':pl.Int64,'cdr3_end':pl.Int64,
 'region_ids':pl.List(pl.Int64),'v_name':pl.String,'j_name':pl.String,'c_name':pl.String,
 'source_contig_id':pl.String,'coding_frame':pl.Int64,'coding_aa_start':pl.Int64,
 'coding_aa_end':pl.Int64,'coverage_kind':pl.String,'sequence_origin':pl.String,
 'source_nt_start':pl.Int64,'source_nt_end':pl.Int64}
_SEARCH_SCHEMA={'entry_id':pl.String,'source_contig_id':pl.String,'locus':pl.String,
 'seq':pl.String,'coding_frame':pl.Int64,'coding_aa_start':pl.Int64,'coding_aa_end':pl.Int64,
 'source_nt_start':pl.Int64,'source_nt_end':pl.Int64}
_CATALOG_SCHEMA={'sequence_id':pl.String,'locus':pl.String,'raw_nt_length':pl.Int64,
 'rev_comp':pl.Boolean,'coding_frame':pl.Int64,'coding_aa_length':pl.Int64,
 'supported_annotated_aa':pl.Int64,'status':pl.String,'support_semantics':pl.String}


def adapt_contigs(raw_records, annotations, supports, *, boundaries=None):
    """Return catalog, heavy candidates and all-locus supported search segments.

    ``supports`` maps sequence_id to sequence and retained_primary_exact_match_mask.
    Native AIRR coordinates are 1-based closed, on the coding strand; returned
    fragment offsets are 0-based half-open amino-acid coordinates. False masks
    mean no retained-primary exact match, not proven biological absence.
    """
    raw_records=list(raw_records);annotations=list(annotations)
    def keyed(rows):
        result={r['sequence_id']:r for r in rows}
        if len(result)!=len(rows):raise ValueError('duplicate sequence_id')
        return result
    raw=keyed(raw_records);native=keyed(annotations)
    projected = keyed(boundaries.to_dicts()) if boundaries is not None else None
    if projected is not None and set(projected) != set(raw):
        raise ValueError('observed boundary keys differ from contigs')
    if set(raw)!=set(native) or set(raw)!=set(supports):
        raise ValueError('raw, annotation and support keys differ')
    catalog=[];candidates=[];entries=[]
    for name,row in raw.items():
        a=native[name];s=supports[name];nt=row['sequence']
        if a.get('sequence')!=nt or s['sequence']!=nt:raise ValueError('sequence source mismatch')
        mask=s['retained_primary_exact_match_mask']
        if len(mask)!=len(nt) or any(type(v) is not bool for v in mask):
            raise ValueError('invalid retained-primary support mask')
        if not nt or set(nt)-set('ACGTN'):raise ValueError('invalid raw nucleotide alphabet')
        rev=a.get('rev_comp')
        locus=a.get('locus') or 'unknown'
        item=dict(sequence_id=name,locus=locus,raw_nt_length=len(nt),
                  rev_comp=None if rev in (None,'') else rev=='T',
                  coding_frame=None,coding_aa_length=None,supported_annotated_aa=0,
                  status='missing_or_ambiguous_junction_frame',
                  support_semantics='retained_primary_exact_matches_not_unique_origin')
        if rev in (None,''):
            item['status']='unavailable_coding_orientation'
            catalog.append(item)
            continue
        if rev not in ('T','F'):raise ValueError('missing coding orientation')
        if rev=='T':nt=reverse_complement(nt);mask=mask[::-1]
        junction=a.get('junction_aa') or ''
        # A native junction uniquely locating the translation is required; do not
        # infer a reading frame from planted sequences or a germline completion.
        options=[]
        if len(junction)>=3 and junction[0]=='C' and junction[-1] in 'FW':
            for frame in range(3):
                aa=translate(nt,frame=frame)
                if aa.count(junction)==1:options.append((frame,aa,aa.index(junction)))
        if len(options)!=1:catalog.append(item);continue
        frame,aa,jstart=options[0];jend=jstart+len(junction)
        labels=[-1]*len(aa)
        for region,label in [('fwr1',0),('cdr1',1),('fwr2',2),('cdr2',3),('fwr3',4),('cdr3',5),('fwr4',7)]:
            start,end=a.get(region+'_start'),a.get(region+'_end')
            if start in (None,'') or end in (None,''):continue
            start,end=int(start),int(end)
            if not 1<=start<=end<=len(nt):raise ValueError('native region outside query')
            lo=max(0,(start-1-frame+2)//3);hi=min(len(aa),(end-frame)//3)
            labels[lo:hi]=[label]*max(0,hi-lo)
        labels[jstart:jend]=[4]+[5]*(len(junction)-2)+[7]
        if projected is not None and projected[name]['status'] == 'resolved':
            boundary = projected[name]['mature_start']
            if type(boundary) is not int or not 1 <= boundary <= len(nt):
                raise ValueError('observed mature boundary outside contig')
            lo=max(0,(boundary-1-frame+2)//3)
            labels[:lo]=[-1]*lo
            start,end=a.get('fwr1_start'),a.get('fwr1_end')
            if start not in (None,'') and end not in (None,''):
                hi=max(0,(int(start)-1-frame+2)//3)
                if boundary <= int(start):labels[lo:hi]=[0]*max(0,hi-lo)
        # Only native alignment-supported constant positions, never an appended reference.
        if a.get('c_call') and a.get('fwr4_end') and float(a.get('mmseqs2_tend') or 0)>float(a.get('mmseqs2_t_vjend') or 0)>0:
            lo=max(0,(int(a['fwr4_end'])-frame)//3)
            hi=min(len(aa),(int(float(a.get('mmseqs2_qend') or 0))-frame)//3)
            for pos in range(lo,hi):
                if labels[pos]<0:labels[pos]=8
        good=[labels[i]>=0 and aa[i] in 'ACDEFGHIKLMNPQRSTVWY' and
              all(mask[frame+3*i:frame+3*i+3]) for i in range(len(aa))]
        runs=[];offset=0
        for eligible,group in groupby(good):
            length=sum(1 for _ in group)
            if eligible and length>=MIN_QUERY_LEN:runs.append((offset,offset+length))
            offset+=length
        for n,(lo,hi) in enumerate(runs):
            entries.append(dict(entry_id=f'{name}_segment{n}',source_contig_id=name,locus=locus,
                                seq=aa[lo:hi],coding_frame=frame,coding_aa_start=lo,coding_aa_end=hi,
                                source_nt_start=len(nt)-(frame+3*hi) if rev=='T' else frame+3*lo,
                                source_nt_end=len(nt)-(frame+3*lo) if rev=='T' else frame+3*hi))
        choices=[(lo,hi) for lo,hi in runs if lo<=jstart and jend<=hi]
        item.update(coding_frame=frame,coding_aa_length=len(aa),supported_annotated_aa=sum(good),
                    status='non_heavy_observation' if locus!='IGH' else 'no_supported_complete_junction')
        if locus=='IGH' and choices:
            lo,hi=max(choices,key=lambda p:(p[1]-p[0],-p[0]))
            candidates.append(dict(clone_id=name,seq=aa[lo:hi],nt_seq=nt[frame+3*lo:frame+3*hi],
              cdr3_seq=junction,cdr3_start=jstart-lo,cdr3_end=jend-lo,region_ids=labels[lo:hi],
              v_name=a.get('v_call') or 'unknown',j_name=a.get('j_call') or 'unknown',
              c_name=a.get('c_call') or 'unknown',source_contig_id=name,coding_frame=frame,
              coding_aa_start=lo,coding_aa_end=hi,coverage_kind='supported_contig_fragment',
              sequence_origin='native_preimputation_contig_retained_primary_support',
              source_nt_start=len(nt)-(frame+3*hi) if rev=='T' else frame+3*lo,
              source_nt_end=len(nt)-(frame+3*lo) if rev=='T' else frame+3*hi))
            item['status']='usable_heavy_candidate'
        catalog.append(item)
    return {'catalog':pl.DataFrame(catalog,schema=_CATALOG_SCHEMA),
            'candidates':pl.DataFrame(candidates,schema=_CANDIDATE_SCHEMA),
            'search_entries':pl.DataFrame(entries,schema=_SEARCH_SCHEMA)}
