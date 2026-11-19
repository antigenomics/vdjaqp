"""Protein-level isotype evidence from peptides that span the J–C junction.

A peptide placement covering at least one FR4 residue (``region_ids`` 7)
immediately followed by a native constant residue (``region_ids`` 8) links the
candidate's V(D)J identity to a constant region in one observed peptide. Its
constant-side residues are compared with the human IGHC CH1 references below; the
isotype is called from those residues and then compared with the candidate's
``c_name``. Nothing here appends reference constant sequence to candidates:
only constant positions that the candidate already carries can be spanned
(see :func:`vdjaqp.contigs.adapt_contigs`).

Coordinates are zero-based half-open placement offsets in the candidate ``seq``.
The first constant residue after FR4 is CH1 residue 1 (the split codon completed
by the J terminal G). Residues are compared under I/L equivalence.
"""
from __future__ import annotations

import polars as pl

#: Human IGHC CH1 amino-acid references (heavy chain only). Translated from
#: ARDA ``database/c_genes/human.fasta`` (CH1 exons lifted from the genome
#: assembly; SHA256 982d81a4...6362) after prefixing the IGHJ terminal G that
#: completes the split first codon. That record lacks the first CH1 base of
#: IGHG4, so its translation starts at residue 2 and residue 1 (Ala) is the
#: completed split codon. IGHM equals ENSP00000375001.2 over all 104 residues;
#: IGHG1 equals UniProt P01857-1 except residue 97 (R here, K in P01857-1).
IGHC_CH1 = {
    'IGHA1': 'ASPTSPKVFPLSLCSTQPDGNVVIACLVQGFFPQEPLSVTWSESGQGVTARNFPPSQDASGDLYTTSSQLTLPATQCLAGKSVTCHVKHYTNPSQDVTVPCP',
    'IGHA2': 'ASPTSPKVFPLSLDSTPQDGNVVVACLVQGFFPQEPLSVTWSESGQNVTARNFPPSQDASGDLYTTSSQLTLPATQCPDGKSVTCHVKHYTNPSQDVTVPC',
    'IGHD': 'APTKAPDVFPIISGCRHPKDNSPVVLACLITGYHPTSVTVTWYMGTQSQPQRTFPEIQRRDSYYMTSSQLSTPLQQWRQGEYKCVVQHTASKSKKEIFRWP',
    'IGHE': 'ASTQSPSVFPLTRCCKNIPSNATSVTLGCLATGYFPEPVMVTWDTGSLNGTTMTLPATTLTLSGHYATISLLTVSGAWAKQMFTCRVAHTPSSTDWVDNKTFS',
    'IGHG1': 'ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKRV',
    'IGHG2': 'ASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVTSSNFGTQTYTCNVDHKPSNTKVDKTV',
    'IGHG3': 'ASTKGPSVFPLAPCSRSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYTCNVNHKPSNTKVDKRV',
    'IGHG4': 'ASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTKTYTCNVDHKPSNTKVDKRV',
    'IGHM': 'GSASAPTLFPLVSCENSPSDTSSVAVGCLAQDFLPDSITFSWKYKNNSDISSTRGFPSVLRGGKYAATSQVLLPSKDVMQGTDEHVVCKVQHPNGNKEKNVPLP',
}
IGHC_CH1_SOURCE = ('ARDA database/c_genes/human.fasta CH1 exons (sha256 '
                   '982d81a46187a4664f31dce37dc61e235192a80e89c2ad92382bc936f0ca6362), '
                   'translated with the IGHJ terminal G completing codon 1')

_FR4, _C = 7, 8
_ISOTYPE = {'IGHG': 'IgG', 'IGHA': 'IgA', 'IGHM': 'IgM', 'IGHD': 'IgD', 'IGHE': 'IgE'}

_PLACEMENT_SCHEMA = {'clone_id': pl.String, 'peptide_id': pl.String, 'start': pl.Int64,
                     'end': pl.Int64, 'c_offset': pl.Int64, 'c_residues': pl.String,
                     'consistent_genes': pl.String, 'consistent_isotypes': pl.String}


def isotype_of_gene(gene: str | None) -> str | None:
    """``IGHG1*01`` -> ``IgG``; unknown, empty or non-IGH calls -> ``None``."""
    if not gene:
        return None
    return _ISOTYPE.get(gene.strip().split('*')[0][:4])


def c_name_isotype(c_name: str | None) -> str | None:
    """Isotype of a (possibly comma-separated) ``c_name``; ``None`` if absent or mixed."""
    if not c_name:
        return None
    calls = {isotype_of_gene(g) for g in c_name.split(',')}
    return calls.pop() if len(calls) == 1 and None not in calls else None


def _il(s: str) -> str:
    return s.replace('I', 'L')


def _class_key(candidates: pl.DataFrame) -> pl.Expr:
    if 'identity_class' in candidates.columns:
        return pl.col('identity_class')
    return pl.col('cdr3_seq').str.replace_all('I', 'L')


def junction_constant_placements(placements: pl.DataFrame, candidates: pl.DataFrame,
                                 references: dict[str, str] | None = None) -> pl.DataFrame:
    """One row per placement whose interval contains an FR4→C boundary.

    ``placements`` needs ``clone_id``, ``peptide_id``, ``start`` and ``end``; an
    optional ``peptide`` column supplies the observed residues (otherwise the
    candidate residues under the placement are used). ``candidates`` needs
    ``clone_id``, ``seq`` and ``region_ids``. ``c_offset`` is the CH1 residue
    index (0-based) of the first constant residue covered. ``consistent_genes``
    lists every reference whose CH1 residues equal the covered constant residues
    under I/L equivalence; an empty string means none agrees.
    """
    references = IGHC_CH1 if references is None else references
    refs = {g: _il(s) for g, s in references.items()}
    with_c = candidates.filter(pl.col('region_ids').list.contains(_C)).select('clone_id', 'seq', 'region_ids')
    if with_c['clone_id'].n_unique() != with_c.height:
        raise ValueError('duplicate candidate clone_id')
    cols = ['clone_id', 'peptide_id', 'start', 'end'] + (['peptide'] if 'peptide' in placements.columns else [])
    hits = placements.select(cols).join(with_c, on='clone_id', how='inner')
    rows = []
    for r in hits.iter_rows(named=True):
        start, end, labels, seq = r['start'], r['end'], r['region_ids'], r['seq']
        if not 0 <= start < end <= len(seq) or len(labels) != len(seq):
            raise ValueError('placement outside candidate sequence: ' + r['clone_id'])
        peptide = r.get('peptide') or seq[start:end]
        if len(peptide) != end - start:
            raise ValueError('peptide length differs from placement interval: ' + r['peptide_id'])
        boundary = next((i for i in range(start, end - 1)
                         if labels[i] == _FR4 and labels[i + 1] == _C), None)
        if boundary is None:
            continue
        c0 = boundary + 1
        c_first = c0
        while c_first > 0 and labels[c_first - 1] == _C:
            c_first -= 1
        offset = c0 - c_first
        stop = c0
        while stop < end and labels[stop] == _C:
            stop += 1
        c_res = _il(peptide[c0 - start:stop - start])
        genes = sorted(g for g, s in refs.items() if s[offset:offset + len(c_res)] == c_res)
        rows.append(dict(clone_id=r['clone_id'], peptide_id=r['peptide_id'], start=start, end=end,
                         c_offset=offset, c_residues=c_res, consistent_genes=','.join(genes),
                         consistent_isotypes=','.join(sorted({isotype_of_gene(g) or g for g in genes}))))
    return pl.DataFrame(rows, schema=_PLACEMENT_SCHEMA)


def _call(consistent: pl.DataFrame) -> tuple[str | None, str | None]:
    specific = consistent.filter(~pl.col('consistent_isotypes').str.contains(','))
    isotypes = set(specific['consistent_isotypes'])
    if len(isotypes) != 1:
        return ('conflicting' if isotypes else ('unresolved' if consistent.height else None)), None
    isotype = isotypes.pop()
    genes = None
    for g in specific.filter(pl.col('consistent_isotypes') == isotype)['consistent_genes']:
        genes = set(g.split(',')) if genes is None else genes & set(g.split(','))
    return isotype, (genes.pop() if genes and len(genes) == 1 else None)


def isotype_evidence(placements: pl.DataFrame, candidates: pl.DataFrame,
                     scores: pl.DataFrame | None = None, *, threshold: float = 20.0,
                     references: dict[str, str] | None = None) -> pl.DataFrame:
    """Per-identity-class protein-level isotype evidence.

    Classes are ``identity_class`` when ``candidates`` carries it, else the I/L
    canonical ``cdr3_seq`` (anchored junction). Every class in ``candidates`` is
    returned, including those without constant residues or J–C placements.

    Columns:
        ``jc_placements`` / ``jc_peptides``: J–C spanning placements and distinct
        peptide IDs on the class's records. ``jc_peptides_unique_to_class``: those
        whose placements (any region, all supplied placements) fall on this class
        only, so they link this V(D)J identity to the constant region. ``reference_consistent_peptides``:
        J–C peptides whose constant residues equal at least one CH1 reference.
        ``protein_isotype``: the single isotype supported by isotype-specific
        consistent peptides; ``'conflicting'`` if several isotypes are each
        specifically supported, ``'unresolved'`` if consistent peptides fit
        several isotypes, null without consistent J–C evidence.
        ``protein_subclass``: the gene when every supporting peptide fits only it,
        else null (an isotype-level call, e.g. a CH1 prefix shared by IGHG1–4).
        ``identity_linked_isotype``: the same call from class-unique J–C peptides only.
        ``c_name_isotype``: isotype of the records' ``c_name`` (null if absent or
        mixed). ``c_name_agreement``: ``agree``/``disagree``/``no_c_name``/``not_called``.
        ``accepted``: any record with ``confirmed_nats >= threshold`` (null without
        ``scores``).
    """
    key = _class_key(candidates)
    records = candidates.select('clone_id', key.alias('identity_class'), 'c_name',
                                pl.col('region_ids').list.contains(_C).alias('has_constant'))
    jc = junction_constant_placements(placements, candidates, references).join(
        records.select('clone_id', 'identity_class'), on='clone_id', how='left')
    # Class specificity counts every placement of the peptide, not only J–C ones.
    spread = placements.select('clone_id', 'peptide_id').join(
        jc.select('peptide_id').unique(), on='peptide_id').join(
        records.select('clone_id', 'identity_class'), on='clone_id', how='left').group_by(
        'peptide_id').agg(pl.col('identity_class').n_unique().alias('n_classes'))
    jc = jc.join(spread, on='peptide_id')
    classes = records.group_by('identity_class', maintain_order=True).agg(
        pl.len().alias('records'), pl.col('has_constant').sum().alias('records_with_constant'),
        pl.col('c_name').fill_null('unknown').unique().sort().str.join(',').alias('c_names'))
    per = {}
    for cls, group in jc.group_by('identity_class'):
        cls = cls[0]
        consistent = group.filter(pl.col('consistent_genes') != '')
        isotype, subclass = _call(consistent)
        linked, _ = _call(consistent.filter(pl.col('n_classes') == 1))
        per[cls] = dict(identity_class=cls, jc_placements=group.height,
                        jc_peptides=group['peptide_id'].n_unique(),
                        jc_peptides_unique_to_class=group.filter(pl.col('n_classes') == 1)['peptide_id'].n_unique(),
                        reference_consistent_peptides=consistent['peptide_id'].n_unique(),
                        protein_isotype=isotype, protein_subclass=subclass, identity_linked_isotype=linked)
    schema = {'identity_class': pl.String, 'jc_placements': pl.Int64, 'jc_peptides': pl.Int64,
              'jc_peptides_unique_to_class': pl.Int64, 'reference_consistent_peptides': pl.Int64,
              'protein_isotype': pl.String, 'protein_subclass': pl.String,
              'identity_linked_isotype': pl.String}
    evidence = pl.DataFrame(list(per.values()), schema=schema)
    out = classes.join(evidence, on='identity_class', how='left').with_columns(
        pl.col('jc_placements', 'jc_peptides', 'jc_peptides_unique_to_class',
               'reference_consistent_peptides').fill_null(0),
        pl.col('c_names').map_elements(c_name_isotype, return_dtype=pl.String).alias('c_name_isotype'))
    out = out.with_columns(
        pl.when(pl.col('protein_isotype').is_null() | pl.col('protein_isotype').is_in(['conflicting', 'unresolved']))
        .then(pl.lit('not_called'))
        .when(pl.col('c_name_isotype').is_null()).then(pl.lit('no_c_name'))
        .when(pl.col('protein_isotype') == pl.col('c_name_isotype')).then(pl.lit('agree'))
        .otherwise(pl.lit('disagree')).alias('c_name_agreement'))
    if scores is None:
        return out.with_columns(pl.lit(None, pl.Float64).alias('max_confirmed_nats'),
                                pl.lit(None, pl.Boolean).alias('accepted'))
    best = scores.select('clone_id', 'confirmed_nats').join(
        records.select('clone_id', 'identity_class'), on='clone_id', how='inner').group_by(
        'identity_class').agg(pl.col('confirmed_nats').max().alias('max_confirmed_nats'))
    return out.join(best, on='identity_class', how='left').with_columns(
        (pl.col('max_confirmed_nats').fill_null(float('-inf')) >= threshold).alias('accepted'))
