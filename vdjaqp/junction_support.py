"""Contextual, inferred-N support from existing placements; no search or scoring."""
import numpy as np
import polars as pl
from arda.refbuild.translate import translate
from .match import ms_collapse


def _keyed(frame, key):
    rows=frame.to_dicts(); result={r[key]:r for r in rows}
    if len(rows)!=len(result) or any(not isinstance(k,str) or not k for k in result):
        raise ValueError('catalogue keys must be unique nonempty strings')
    return result


def observed_junction_classes(candidates):
    """Complete concordant anchored observed junctions; other candidates stay unresolved."""
    result={}
    for key,c in _keyed(candidates,'clone_id').items():
        junction,offset,seq=c.get('cdr3_seq'),c.get('cdr3_start'),c['seq']
        if (isinstance(junction,str) and len(junction)>=2 and junction[0]=='C' and
            junction[-1] in 'FW' and not set(junction)-set('ACDEFGHIKLMNPQRSTVWY') and
            type(offset) is int and offset>=0 and seq[offset:offset+len(junction)]==junction):
            result[key]=ms_collapse([junction])[0]
    return result


def inferred_n_footprints(candidates, annotations):
    """Return per-candidate masks from verified native inferred nucleotide geometry.

    Candidate columns: clone_id, seq, cdr3_start (0-based AA offset), cdr3_seq
    (anchored junction AA). Annotation columns: clone_id, cdr3_nt, v_end_nt,
    j_start_nt, d_start_nt, d_end_nt, np1, np2. Native N slices are exactly
    nt[v_end_nt:d_start_nt-1] and nt[d_end_nt:j_start_nt]; D is 1-based closed.
    Missing annotations/coordinates remain unavailable. Never fills absent geometry.
    """
    catalog=_keyed(candidates,'clone_id'); annotation=_keyed(annotations,'clone_id')
    if set(annotation)-set(catalog):raise ValueError('annotation has unknown candidate')
    result={};observed=observed_junction_classes(candidates)
    for key,c in catalog.items():
        seq,junction,offset=c['seq'],c.get('cdr3_seq'),c.get('cdr3_start')
        if not isinstance(seq,str) or not seq or set(seq)-set('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError('invalid candidate sequence')
        masks={name:np.zeros(len(seq),dtype=bool) for name in ('n','mixed','d')}
        a=annotation.get(key);status='missing_geometry'
        if a is not None:
            fields=[a.get(k) for k in ('v_end_nt','j_start_nt','d_start_nt','d_end_nt')]
            nt=a.get('cdr3_nt')
            if nt is not None and all(x is not None for x in fields):
                if key not in observed:raise ValueError('inferred geometry needs complete concordant observed candidate junction')
                if not isinstance(nt,str) or not nt or set(nt)-set('ACGT') or any(type(x) is not int for x in fields):
                    raise ValueError('invalid inferred nucleotide geometry')
                lo,hi,ds,de=fields
                if not 0<=lo<=ds-1<de<=hi<=len(nt):raise ValueError('inferred tract bounds overlap or leave junction')
                if (len(nt)!=3*len(junction) or translate(nt)!=junction or
                    (a.get('np1') or '')!=nt[lo:ds-1] or (a.get('np2') or '')!=nt[de:hi]):
                    raise ValueError('inferred sequence, N slices or translation disagree')
                starts=3*np.arange(len(junction));ends=starts+3
                whole=np.zeros(len(junction),dtype=bool);overlap=whole.copy()
                for left,right in ((lo,ds-1),(de,hi)):
                    whole|=(starts>=left)&(ends<=right)
                    overlap|=(starts<right)&(ends>left)
                masks['n'][offset:offset+len(junction)]=whole
                masks['mixed'][offset:offset+len(junction)]=overlap&~whole
                masks['d'][offset:offset+len(junction)]=(starts<de)&(ends>ds-1)
                status='inferred_geometry'
        result[key]=dict(status=status,**masks)
    return result


def contextual_n_support(candidates, peptides, placements, annotations, *,
                         accepted_scores, accepted_support, min_consecutive=3):
    """Keep every placement, labelling rescue eligibility without changing scores.

    Peptides: peptide_id, seq. Placements: observation_id, peptide_id, clone_id,
    start, end (0-based half-open candidate AA), is_germline (bool).
    accepted_support: observation_id, clone_id, explicitly credited support from
    the accepted scorer (not a list of arbitrary compatible placements).
    accepted_scores: clone_id, confirmed_nats from that same accepted output.
    """
    if type(min_consecutive) is not int or min_consecutive<1:
        raise ValueError('min_consecutive must be a positive integer')
    cs=_keyed(candidates,'clone_id');ps=_keyed(peptides,'peptide_id')
    footprints=inferred_n_footprints(candidates,annotations)
    scores=_keyed(accepted_scores,'clone_id')
    if set(scores)-set(cs):raise ValueError('score has unknown candidate')
    for r in scores.values():
        if isinstance(r['confirmed_nats'],bool) or not np.isfinite(r['confirmed_nats']):raise ValueError('invalid accepted score')
    explained=set()
    for r in accepted_support.to_dicts():
        if r['clone_id'] not in scores or not isinstance(r['observation_id'],str) or not r['observation_id']:
            raise ValueError('credited support requires known score and observation')
        if scores[r['clone_id']]['confirmed_nats']>=20:explained.add(r['observation_id'])
    encoded={k:np.frombuffer(ms_collapse([r['seq']])[0].encode(),dtype=np.uint8) for k,r in cs.items()}
    queries={}
    for k,r in ps.items():
        if not isinstance(r['seq'],str) or not r['seq'] or set(r['seq'])-set('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError('invalid peptide sequence')
        queries[k]=np.frombuffer(ms_collapse([r['seq']])[0].encode(),dtype=np.uint8)
    rows=[]
    # ponytail: per-placement dispatch; use chunked array reductions if profiling identifies a bottleneck.
    for index,r in enumerate(placements.to_dicts()):
        c,p,o=r['clone_id'],r['peptide_id'],r['observation_id'];start,end=r['start'],r['end']
        if c not in cs or p not in ps or not isinstance(o,str) or not o or type(r['is_germline']) is not bool:
            raise ValueError('invalid placement keys or germline label')
        if type(start) is not int or type(end) is not int or not 0<=start<end<=len(encoded[c]) or end-start!=len(queries[p]):
            raise ValueError('invalid whole-peptide placement')
        f=footprints[c];agree=encoded[c][start:end]==queries[p];support=agree&f['n'][start:end]
        edges=np.flatnonzero(np.diff(np.r_[False,support,False].astype(np.int8)))
        longest=int((edges[1::2]-edges[::2]).max(initial=0))
        available=f['status']=='inferred_geometry';excluded=o in explained
        rows.append(dict(placement_index=index,observation_id=o,peptide_id=p,clone_id=c,
            start=start,end=end,is_germline=r['is_germline'],geometry_status=f['status'],
            geometry_provenance='inferred nucleotide reconstruction; not measured N bases',
            minimum_run=min_consecutive,longest_n_run=longest,
            agreeing_n_residues=int(support.sum()),mixed_codon_residues=int(f['mixed'][start:end].sum()),
            d_overlapping_residues=int(f['d'][start:end].sum()),excluded_mapped_observation=excluded,
            rescue_eligible=available and not excluded and not r['is_germline'] and longest>=min_consecutive))
    return pl.DataFrame(rows) if rows else pl.DataFrame(schema={'placement_index':pl.Int64,
        'observation_id':pl.String,'peptide_id':pl.String,'clone_id':pl.String,
        'start':pl.Int64,'end':pl.Int64,'is_germline':pl.Boolean,'geometry_status':pl.String,
        'geometry_provenance':pl.String,'minimum_run':pl.Int64,'longest_n_run':pl.Int64,
        'agreeing_n_residues':pl.Int64,'mixed_codon_residues':pl.Int64,
        'd_overlapping_residues':pl.Int64,'excluded_mapped_observation':pl.Boolean,'rescue_eligible':pl.Boolean})
