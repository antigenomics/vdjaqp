"""Explicit known-panel insertion for experiment-matched reconstruction only."""
import polars as pl
from arda.refbuild.translate import translate
from .adapter import to_airr_compatible
from .airr_ms_inputs import attach_background_membership


def append_known_references(candidates, search_entries, reference_panel, *,
                            expected_count, source_format):
    """Append supplied annotated heavy references without claiming RNA support.

    Background tables are unchanged native ARDA FASTQ exports. The supplied
    panel uses scorer columns, including coding ``nt_seq``; no annotation, truth
    selection, germline completion or sequence generation happens here. Identical
    sequences with different IDs remain distinct candidate rows intentionally.
    """
    if type(expected_count) is not int or expected_count < 1 or reference_panel.height != expected_count:
        raise ValueError('reference count must equal the declared positive panel count')
    candidates=attach_background_membership(candidates,inserted_reference_ids=[])
    panel = to_airr_compatible(reference_panel, source_format=source_format)
    to_airr_compatible(candidates, source_format='arda_fastq_native')
    if 'nt_seq' not in reference_panel.columns:
        raise ValueError('known references require their supplied coding nucleotide sequence')
    if 'candidate_construction' in candidates.columns or ('sequence_origin' in search_entries.columns and search_entries.filter(pl.col('sequence_origin')=='supplied_known_reference').height):
        raise ValueError('background must be the unmodified observation-only adapter output')
    required={'entry_id','seq','locus','source_observation_id','parent_clone_id'}
    if not required <= set(search_entries.columns) or any(search_entries.schema[k]!=pl.String for k in required):
        raise ValueError('search entries require native string sequence/source fields')
    if search_entries.select((~pl.col('seq').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]{6,}$') |
                              pl.col('source_observation_id').is_null() | pl.col('parent_clone_id').is_null() |
                              ~pl.col('locus').is_in(['IGH','IGK','IGL'])).fill_null(True).any()).item():
        raise ValueError('search entries require canonical observed sequences and native observation IDs')
    origins=['arda_observed_read','arda_native_assembled_read','arda_corrected_junction']
    required_native={'sequence_origin','coverage_kind','nt_seq','native_junction_nt',
                     'native_junction_aa','native_clone_row','parent_clone_id'}
    if not required_native <= set(candidates.columns):
        raise ValueError('background candidates require native ARDA sequence provenance')
    for row in candidates.iter_rows(named=True):
        origin=row['sequence_origin']
        if (origin not in origins or set(row['seq'])-set('ACDEFGHIKLMNPQRSTVWY') or
            not row['cdr3_seq'].startswith('C') or row['cdr3_seq'][-1] not in 'FW' or
            row['coverage_kind'] != (
                'junction_only' if origin=='arda_corrected_junction' else 'observed_variable_fragment') or
            not isinstance(row['nt_seq'],str) or set(row['nt_seq'])-set('ACGT') or
            translate(row['nt_seq'])!=row['seq'] or len(row['nt_seq'])!=3*len(row['seq']) or
            row['native_junction_aa']!=row['cdr3_seq'] or
            not isinstance(row['native_junction_nt'],str) or
            set(row['native_junction_nt'])-set('ACGT') or
            len(row['native_junction_nt'])!=3*len(row['native_junction_aa']) or
            translate(row['native_junction_nt'])!=row['native_junction_aa'] or
            type(row['native_clone_row']) is not int or row['native_clone_row']<0 or
            not isinstance(row['parent_clone_id'],str) or not row['parent_clone_id'] or
            (origin=='arda_corrected_junction' and row['seq']!=row['native_junction_aa'])):
            raise ValueError('background candidates must retain valid native ARDA coding/provenance fields')
    for key,frame in [('clone_id',candidates),('entry_id',search_entries)]:
        if key not in frame.columns or frame[key].null_count() or frame[key].n_unique()!=frame.height:
            raise ValueError('background identifiers must be unique and non-null')
    occupied=set(candidates['clone_id']) | set(search_entries['entry_id'])
    if occupied & set(reference_panel['clone_id']):
        raise ValueError('reference/background identifier collision')
    for row in reference_panel.iter_rows(named=True):
        nt,aa,junction=row['nt_seq'],row['seq'],row['cdr3_seq']
        if (not isinstance(nt,str) or set(nt)-set('ACGT') or
            not 3*len(aa)<=len(nt)<3*len(aa)+3 or translate(nt)!=aa or
            set(aa)-set('ACDEFGHIKLMNPQRSTVWY') or junction[0]!='C' or junction[-1] not in 'FW' or
            not (row['v_name'] or '').startswith('IGHV') or not (row['j_name'] or '').startswith('IGHJ')):
            raise ValueError('known panel must contain valid coding heavy sequences and anchored junctions')
    fields=['clone_id','seq','nt_seq','cdr3_seq','cdr3_start','cdr3_end',
            'region_ids','v_name','j_name','c_name']
    inserted=reference_panel.select(fields).with_columns(
        pl.lit(False).alias('background_member'),
        pl.lit('reference_panel_insertion').alias('candidate_construction'),
        pl.lit('supplied_known_reference').alias('sequence_origin'),
        pl.lit('supplied_reference_sequence_not_rna_support').alias('coverage_kind'),
        pl.lit(source_format).alias('reference_source_format'))
    added_entries=reference_panel.select(pl.col('clone_id').alias('entry_id'),'seq').with_columns(
        pl.lit('IGH').alias('locus'),
        pl.lit(False).alias('background_member'),
        pl.lit('supplied_known_reference').alias('sequence_origin'))
    return dict(
        candidates=pl.concat([candidates.with_columns(pl.lit('fastq_derived').alias('candidate_construction')),inserted],how='diagonal_relaxed'),
        search_entries=pl.concat([search_entries,added_entries],how='diagonal_relaxed'),
        reference_panel=panel)


if __name__ == '__main__':
    import argparse
    import hashlib
    import json
    from pathlib import Path
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidates',type=Path,required=True)
    p.add_argument('--search-entries',type=Path,required=True)
    p.add_argument('--reference-panel',type=Path,required=True)
    p.add_argument('--expected-count',type=int,default=12)
    p.add_argument('--source-format',required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():p.error('output already exists; choose a fresh directory')
    paths={name:getattr(args,name) for name in ('candidates','search_entries','reference_panel')}
    frames={name:pl.read_parquet(path) for name,path in paths.items()}
    result=append_known_references(**frames,expected_count=args.expected_count,source_format=args.source_format)
    args.output.mkdir(parents=True)
    for name,frame in result.items():frame.write_parquet(args.output/(name+'.parquet'),compression='zstd')
    manifest=dict(construction='known-reference insertion; not a generic FASTQ observation',
        source_format=args.source_format,expected_reference_count=args.expected_count,
        input_rows={name:f.height for name,f in frames.items()},output_rows={name:f.height for name,f in result.items()},
        input_sha256={name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in paths.items()},
        output_sha256={name:hashlib.sha256((args.output/(name+'.parquet')).read_bytes()).hexdigest() for name in result})
    (args.output/'insertion_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
