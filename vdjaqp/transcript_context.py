"""Explicit donor-leader imputation; native mature sequences are never completed."""
from collections import defaultdict
import hashlib
import math


def select_contexts(native, sources, alignments=()):
    """Select strict, conservative group, then aligned same-locus donor leaders.

    Native rows: locus, native_allele, gene, core_nt, exact_sources (source IDs).
    Sources: source_id, locus, core_nt, status and validated prefix_nt/leader_nt.
    Alignments: native_allele, source_id, bits. Only validated sources reached by
    strict native mappings are alignment donors. No score threshold is introduced.
    Return flat selections and flat candidate provenance. Missing hits stay explicit.
    """
    native, sources = list(native), list(sources)
    for rows, fields in ((native, ('native_allele', 'gene', 'core_nt', 'locus')),
                         (sources, ('source_id', 'core_nt', 'locus', 'status'))):
        for row in rows:
            if (not isinstance(row, dict) or any(not isinstance(row.get(f), str) or not row[f] for f in fields)
                    or row['locus'] not in ('IGH', 'IGK', 'IGL') or set(row['core_nt']) - set('ACGT')):
                raise ValueError('invalid native/source schema, locus or core')
    for row in native:
        ids = row.get('exact_sources')
        if not isinstance(ids, list) or any(not isinstance(s, str) for s in ids):
            raise ValueError('invalid exact-source identities')
    for row in sources:
        if row['status'] not in ('valid', 'coding_source_coordinate_mismatch',
                'missing_leader_coordinates', 'invalid_leader_coordinates',
                'invalid_leader_frame_start_or_stop', 'leader_core_translation_mismatch',
                'invalid_utr_coordinates', 'invalid_utr_sequence'):
            raise ValueError('unknown source-context status')
    ns = {r['native_allele']: r for r in native}
    ss = {r['source_id']: r for r in sources}
    if len(ns) != len(native) or len(ss) != len(sources):
        raise ValueError('duplicate native allele or source ID')
    groups = defaultdict(set)
    strict = {}
    for source in sources:
        if source['status'] == 'valid':
            from arda.refbuild.translate import translate
            prefix, leader = source.get('prefix_nt'), source.get('leader_nt')
            if (not isinstance(prefix, str) or not isinstance(leader, str) or not leader or
                    not prefix.endswith(leader) or set(prefix + source['core_nt']) - set('ACGT') or
                    not leader.startswith('ATG') or len(leader) % 3 or '*' in translate(leader)):
                raise ValueError('invalid source leader')
    for row in native:
        ids = set(row['exact_sources'])
        for sid in ids:
            if sid not in ss:
                raise ValueError('unknown exact-source identity')
            source = ss[sid]
            if source['locus'] != row['locus'] or source['core_nt'] != row['core_nt']:
                raise ValueError('exact source bridge disagrees with native core/locus')
        groups[row['locus'], row['gene']].update(ids)
        if ids and all(ss[s]['status'] == 'valid' for s in ids) and len({ss[s]['prefix_nt'] for s in ids}) == 1:
            strict[row['native_allele']] = ids
    donors = set().union(*strict.values()) if strict else set()
    hits = defaultdict(list)
    for hit in alignments:
        if not isinstance(hit, dict):
            raise ValueError('invalid alignment schema')
        allele, sid, bits = hit.get('native_allele'), hit.get('source_id'), hit.get('bits')
        if allele not in ns or sid not in donors or ns[allele]['locus'] != ss[sid]['locus']:
            raise ValueError('unknown, non-donor or cross-locus alignment')
        if isinstance(bits, bool) or not isinstance(bits, (int, float)) or not math.isfinite(bits):
            raise ValueError('invalid alignment bit score')
        hits[allele].append((sid, bits))
    result, provenance = [], []
    for row in native:
        allele = row['native_allele']
        group = groups[row['locus'], row['gene']]
        candidates, method, confidence, score = set(), 'UNSUPPORTED', 'no_donor_alignment', None
        if allele in strict:
            candidates, method, confidence = strict[allele], 'STRICT_EXACT', 'validated_unique_context'
        elif group and all(ss[s]['status'] == 'valid' for s in group) and len({ss[s]['prefix_nt'] for s in group}) == 1:
            candidates, method, confidence = group, 'IMPUTED_GROUP', 'group_consensus_assumption'
        elif hits[allele]:
            score = max(bits for _, bits in hits[allele])
            candidates = {sid for sid, bits in hits[allele] if bits == score}
            method = 'IMPUTED_ALIGNMENT'
            confidence = ('ambiguous_tied_prefixes' if len({ss[s]['prefix_nt'] for s in candidates}) > 1
                          else 'nearest_core_assumption')
        chosen = min(candidates, key=lambda s: (hashlib.sha256(ss[s]['prefix_nt'].encode()).hexdigest(), s)) if candidates else None
        source = ss[chosen] if chosen is not None else None
        result.append(dict(locus=row['locus'], native_allele=allele, method=method,
            confidence=confidence, selected_source_id=chosen,
            prefix_nt=source['prefix_nt'] if source else None,
            leader_nt=source['leader_nt'] if source else None,
            prefix_sha256=hashlib.sha256(source['prefix_nt'].encode()).hexdigest() if source else None,
            source_candidate_count=len(candidates), alignment_candidate_count=len(hits[allele]),
            selected_bits=score))
        for sid in sorted(candidates):
            provenance.append(dict(native_allele=allele, source_id=sid, method=method,
                                   selected=sid == chosen, prefix_nt=ss[sid]['prefix_nt']))
    return result, provenance


def leader_rna_template(variable_nt, constant_nt, protein_aa, context, *, v_trim_5, translate_nt=None):
    """Prefix an existing mature substrate without changing or completing a base."""
    if translate_nt is None:
        from arda.refbuild.translate import translate as translate_nt
    if type(v_trim_5) is not int or v_trim_5 != 0:
        raise ValueError('unsupported physical join: native V5 trimming must not be filled')
    leader, prefix = context.get('leader_nt'), context.get('prefix_nt')
    if (not variable_nt or not constant_nt or not leader or not prefix or
            set(variable_nt + constant_nt + prefix) - set('ACGT') or
            not prefix.endswith(leader) or not leader.startswith('ATG') or
            len(leader) % 3 or '*' in translate_nt(leader)):
        raise ValueError('invalid coding substrate or leader context')
    mature = variable_nt + constant_nt
    if len(mature) % 3 or translate_nt(mature) != protein_aa + '*':
        raise ValueError('mature coding substrate does not reproduce protein plus stop')
    utr_length = len(prefix) - len(leader)
    template = prefix + mature
    if translate_nt(template[utr_length:]) != translate_nt(leader) + protein_aa + '*':
        raise ValueError('leader/mature translation mismatch')
    return dict(rnaseq_nt=template, mature_nt_start=len(prefix) + 1,
                leader_nt_start=utr_length + 1, mature_nt_length=len(mature))
