"""Project public mature-start markers through observed coding-strand alignments."""
from arda import _markup
import polars as pl


def _project_best_boundaries(queries, contexts, hits):
    """Native ARDA max whole-target bits ranking, retaining every exact-max tie."""
    import math
    if any(type(h.get('bits')) not in (float, int) or not math.isfinite(h['bits']) for h in hits):
        raise ValueError('invalid boundary alignment bit score')
    if not hits:
        return project_transcript_boundaries(queries, contexts, []).with_columns(
            pl.lit(0, dtype=pl.Int64).alias('returned_alignment_count'),
            pl.lit(0, dtype=pl.Int64).alias('lower_score_alignment_count'),
            pl.lit(None, dtype=pl.Float64).alias('best_bits'))
    frame = pl.DataFrame(hits).with_columns(pl.col('bits').max().over('sequence_id').alias('best_bits'))
    retained = frame.filter(pl.col('bits') == pl.col('best_bits'))
    counts = frame.group_by('sequence_id').agg(
        pl.len().cast(pl.Int64).alias('returned_alignment_count'),
        (pl.col('bits') < pl.col('best_bits')).sum().cast(pl.Int64).alias('lower_score_alignment_count'),
        pl.col('best_bits').first())
    result = project_transcript_boundaries(queries, contexts, retained.to_dicts())
    return result.join(counts, on='sequence_id', how='left', maintain_order='left').with_columns(
        pl.col('returned_alignment_count', 'lower_score_alignment_count').fill_null(0))


def search_transcript_boundaries(raw_records, annotations, supports, contexts, *, threads, alignment_output=None):
    """One native batch against the entire public/model catalogue, never a genotype.

    Coding orientation comes solely from observed ARDA annotation. No simulation
    truth, generated-chain identifiers or mature offsets are accepted here.
    """
    import csv
    from pathlib import Path
    from tempfile import TemporaryDirectory
    import os
    from arda import mmseqs
    from arda.refbuild.translate import reverse_complement
    from .resources import resolve_cpu_budget
    threads = resolve_cpu_budget(threads)
    raw_records, annotations, contexts = list(raw_records), list(annotations), list(contexts)
    ann = {r['sequence_id']: r for r in annotations}
    names = [r['sequence_id'] for r in raw_records]
    if len(set(names)) != len(names) or len(ann) != len(annotations) or set(names) != set(ann) or set(names) != set(supports):
        raise ValueError('raw, annotation and support keys must be unique and equal')
    queries, searchable = [], []
    for row in raw_records:
        name, seq = row['sequence_id'], row['sequence']
        if ann[name].get('sequence') != seq or supports[name].get('sequence') != seq:
            raise ValueError('sequence source mismatch')
        mask = supports[name].get('retained_primary_exact_match_mask')
        rev = ann[name].get('rev_comp')
        if rev not in ('F', 'T', None, ''):
            raise ValueError('invalid observed orientation')
        if rev == 'T':
            seq = reverse_complement(seq)
            mask = mask[::-1] if mask is not None else None
        item = dict(sequence_id=name, sequence=seq, retained_primary_exact_match_mask=mask)
        queries.append(item)
        if rev in ('F', 'T'):
            searchable.append(item)
    # Validate even empty/no-hit inputs before native search.
    project_transcript_boundaries(queries, contexts, [])
    if not searchable:
        return _project_best_boundaries(queries, contexts, [])
    if os.environ.get('ARDA_MMSEQS_SEARCH_OPTS', '').strip():
        raise ValueError('unsupported MMseqs override')
    with TemporaryDirectory(prefix='igseq_boundary_') as tmp:
        tmp = Path(tmp)
        qfile, tfile = tmp/'query.fasta', tmp/'target.fasta'
        qfile.write_text(''.join(f">{r['sequence_id']}\n{r['sequence']}\n" for r in searchable))
        tfile.write_text(''.join(f">{r['target']}\n{r['sequence']}\n" for r in contexts))
        qdb = mmseqs.createdb(qfile, tmp/'query', dbtype=2)
        tdb = mmseqs.createdb(tfile, tmp/'target', dbtype=2)
        result = mmseqs.search(qdb, tdb, tmp/'result', tmp/'work', search_type=3,
            sensitivity=7.0, evalue=1e-3, max_seqs=300, threads=threads, kmer=13, extra=['--strand', '1'])
        table = mmseqs.convertalis(qdb, tdb, result, tmp/'hits.tsv', search_type=3, threads=threads)
        if alignment_output is not None:
            import gzip
            import shutil
            with table.open('rb') as source, gzip.open(alignment_output, 'xb') as destination:
                shutil.copyfileobj(source, destination)
        hits = []
        with table.open() as stream:
            for hit in csv.DictReader(stream, fieldnames=mmseqs.DEFAULT_FORMAT_OUTPUT.split(','), delimiter='\t'):
                hits.append(dict(sequence_id=hit['query'], target=hit['target'], qaln=hit['qaln'], taln=hit['taln'], bits=float(hit['bits']),
                    **{key: int(hit[key]) for key in ('qstart','qend','tstart','tend','qlen','tlen')}))
    return _project_best_boundaries(queries, contexts, hits)


def project_transcript_boundaries(queries, contexts, alignments):
    """Return one boundary/status row for every query, including queries without hits.

    Queries contain ``sequence_id, sequence, retained_primary_exact_match_mask``;
    contexts contain ``target, sequence,
    mature_start``. Alignment rows contain ``sequence_id, target, qstart, qend,
    tstart, tend, qlen, tlen, qaln, taln``. All positions are integer 1-based closed
    nucleotide coordinates. Sequences must already be coding-oriented; reverse
    complementation and raw-source offsets belong to the caller. Public context
    start coordinates and any ARDA frame-trim offsets must be verified upstream.

    All supplied alignments are retained alternatives, with no score filtering.
    Both bases adjacent to the boundary must map to adjacent query bases, and
    every alternative must support the same boundary. Both flanking marker bases
    require retained-primary exact-match support in the same coding orientation.
    Missing masks yield unresolved boundaries. This does not certify translation
    frame, cleavage biology or exhaustive search coverage.
    """
    def keyed(rows, key):
        result = {}
        for row in rows:
            name = row.get(key)
            seq = row.get('sequence')
            if not isinstance(name, str) or not name or name in result:
                raise ValueError(f'invalid or duplicate {key}')
            if not isinstance(seq, str) or not seq or set(seq) - set('ACGTN'):
                raise ValueError('expected ungapped uppercase nucleotide sequence')
            result[name] = row
        return result

    queries = keyed(queries, 'sequence_id')
    contexts = keyed(contexts, 'target')
    for query in queries.values():
        mask = query.get('retained_primary_exact_match_mask')
        if mask is not None and (len(mask) != len(query['sequence']) or
                                 any(type(v) is not bool for v in mask)):
            raise ValueError('invalid retained-primary support mask')
    for context in contexts.values():
        if set(context['sequence']) - set('ACGT'):
            raise ValueError('public context must be unambiguous ACGT')
        p = context.get('mature_start')
        if type(p) is not int or not 2 <= p <= len(context['sequence']):
            raise ValueError('mature_start must have a preceding context base')
    evidence = {name: [] for name in queries}
    for hit in alignments:
        name, target = hit.get('sequence_id'), hit.get('target')
        if name not in queries or target not in contexts:
            raise ValueError('unknown query or context key')
        query, context = queries[name]['sequence'], contexts[target]['sequence']
        for prefix, sequence in [('q', query), ('t', context)]:
            start, end, length = (hit.get(prefix + field) for field in ('start', 'end', 'len'))
            if (any(type(v) is not int for v in (start, end, length)) or
                    length != len(sequence) or not 1 <= start <= end <= length):
                raise ValueError('invalid coding-oriented alignment coordinates')
        qaln, taln = hit.get('qaln'), hit.get('taln')
        if (not isinstance(qaln, str) or not isinstance(taln, str) or
                not qaln or len(qaln) != len(taln) or
                set(qaln + taln) - set('ACGTN-') or
                any(q == t == '-' for q, t in zip(qaln, taln))):
            raise ValueError('invalid aligned nucleotide strings')
        if (qaln.replace('-', '') != query[hit['qstart'] - 1:hit['qend']] or
                taln.replace('-', '') != context[hit['tstart'] - 1:hit['tend']]):
            raise ValueError('alignment disagrees with query or context sequence')
        p = contexts[target]['mature_start']
        left, right = _markup.transfer_regions(
            qaln, taln, hit['qstart'], hit['tstart'], [p - 1, p], [p - 1, p])
        brackets = left[0] > 0 and right[0] == left[0] + 1
        mask = queries[name].get('retained_primary_exact_match_mask')
        boundary = (right[0] if brackets and mask is not None and
                    query[left[0] - 1] in 'ACGT' and query[right[0] - 1] in 'ACGT' and
                    mask[left[0] - 1] and mask[right[0] - 1] else None)
        evidence[name].append((target, boundary, brackets))
    rows = []
    for name, hits in evidence.items():
        boundaries = {boundary for _, boundary, _ in hits if boundary is not None}
        supported = sum(boundary is not None for _, boundary, _ in hits)
        bracketing = sum(brackets for _, _, brackets in hits)
        status = ('no_alignment' if not hits else
                  'conflicting_boundary' if len(boundaries) > 1 else
                  'unsupported_boundary' if supported != len(hits) else 'resolved')
        rows.append(dict(sequence_id=name,
                         mature_start=next(iter(boundaries)) if status == 'resolved' else None,
                         status=status, alignment_count=len(hits),
                         candidate_count=len({target for target, _, _ in hits}),
                         bracketing_alignment_count=bracketing,
                         nonbracketing_alignment_count=len(hits) - bracketing,
                         supported_alignment_count=supported))
    return pl.DataFrame(rows, schema={'sequence_id': pl.String, 'mature_start': pl.Int64,
        'status': pl.String, 'alignment_count': pl.Int64, 'candidate_count': pl.Int64,
        'bracketing_alignment_count': pl.Int64, 'nonbracketing_alignment_count': pl.Int64,
        'supported_alignment_count': pl.Int64})
