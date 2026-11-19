"""Exact, reference-backed IGH V annotations; native generation calls stay intact.

Rebuild from pinned installed raw references per invocation, without alignment.
"""
from contextlib import closing
import csv
import hashlib
from importlib.resources import files
import json
from pathlib import Path
import sqlite3
import zipfile


def exact_v_map(legacy, canonical):
    """Join (allele, gene, ungapped nucleotide) rows in one SQLite batch."""
    with closing(sqlite3.connect(':memory:')) as db:
        for table, rows in [('legacy', legacy), ('canonical', canonical)]:
            db.execute(f'create table {table}(name text primary key, gene text, seq text)')
            db.executemany(f'insert into {table} values (?,?,?)',
                           [(name, gene, seq.upper()) for name, gene, seq in rows])
        matches = db.execute('select l.name,l.gene,l.seq,c.name,c.gene from legacy l '
                             'left join canonical c on l.seq=c.seq order by l.name,c.name').fetchall()
    result = {}
    for name, gene, seq, canonical_name, canonical_gene in matches:
        row = result.setdefault(name, dict(legacy_gene=gene,
            sequence_sha256=hashlib.sha256(seq.encode('ascii')).hexdigest(),
            canonical_alleles=[], canonical_genes=[], canonical_matches=[]))
        if canonical_name is not None:
            row['canonical_alleles'].append(canonical_name)
            row['canonical_matches'].append(dict(allele=canonical_name, gene=canonical_gene))
            row['canonical_genes'].append(canonical_gene)
    for row in result.values():
        row['canonical_genes'] = sorted(set(row['canonical_genes']))
        count = len(row['canonical_genes'])
        row.update(canonical_v_gene=row['canonical_genes'][0] if count == 1 else None,
                   canonical_v_gene_count=count,
                   canonical_v_match_status=('unique_gene' if count == 1 else
                                             'ambiguous_genes' if count else 'no_exact_match'))
    return result


def load_reference_gene_mapping(identity_path=None):
    """Verify raw bytes before loading installed GenAIRR reference metadata."""
    identity_path = identity_path or Path(__file__).resolve().parents[1]/'envs/generation/identity.json'
    identity = json.loads(Path(identity_path).read_text())
    base = files('GenAIRR.data')
    legacy_name = 'HUMAN_IGH_OGRDB.pkl'
    canonical_name = 'human_igh_ogrdb_10_observed_v1.zip'
    legacy_bytes = base.joinpath('builtin_dataconfigs', legacy_name).read_bytes()
    canonical_bytes = base.joinpath('published_references', canonical_name).read_bytes()
    hashes = {legacy_name: hashlib.sha256(legacy_bytes).hexdigest(),
              canonical_name: hashlib.sha256(canonical_bytes).hexdigest()}
    expected = {legacy_name: identity['genairr_reference_sha256'][legacy_name],
                canonical_name: identity['canonical_gene_reference_sha256'][canonical_name]}
    for name, digest in hashes.items():
        if digest != expected[name]:
            raise ValueError(f'gene mapping raw reference hash mismatch: {name}')
    # Only load this trusted installed pickle after its pinned digest is checked.
    import pickle
    import io
    cfg = pickle.loads(legacy_bytes)
    with zipfile.ZipFile(io.BytesIO(canonical_bytes)) as archive:
        canonical = json.loads(archive.read('reference.genairr-ref.json'))['catalogue']['alleles']
        source_lock = json.loads(archive.read('source.lock.json'))
    legacy = [(a.name, a.gene, a.ungapped_seq) for group in cfg.v_alleles.values() for a in group]
    canonical = [(a['label'], a['gene'], a['sequence']) for a in canonical if a['segment'] == 'V']
    return dict(method='exact uppercase ungapped nucleotide equality; no alignment or group renaming',
                locus='IGH', segment='V', reference_sha256=hashes,
                canonical_source_lock=source_lock, legacy_allele_count=len(legacy),
                canonical_allele_count=len(canonical), alleles=exact_v_map(legacy, canonical))


def canonical_v_fields(locus, v_call, mapping):
    """Three scalar chain fields; unavailable loci are explicitly unevaluated."""
    if locus != 'IGH':
        return dict(canonical_v_gene=None, canonical_v_match_status='not_evaluated_locus',
                    canonical_v_gene_count=None)
    row = mapping['alleles'][v_call]  # Unknown native calls must fail rather than masquerade as nonmatches.
    return {key: row[key] for key in ('canonical_v_gene', 'canonical_v_match_status', 'canonical_v_gene_count')}


def write_reference_gene_mapping(path, mapping):
    """One exact allele/gene pair per row; explicit nulls for unmatched alleles."""
    fields = ['legacy_allele', 'legacy_gene', 'sequence_sha256', 'canonical_allele',
              'canonical_gene', 'canonical_v_match_status', 'canonical_v_gene_count']
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for name, row in mapping['alleles'].items():
            for match in row['canonical_matches'] or [dict(allele=None, gene=None)]:
                writer.writerow(dict(legacy_allele=name, legacy_gene=row['legacy_gene'],
                    sequence_sha256=row['sequence_sha256'], canonical_allele=match['allele'] or 'null',
                    canonical_gene=match['gene'] or 'null',
                    canonical_v_match_status=row['canonical_v_match_status'],
                    canonical_v_gene_count=row['canonical_v_gene_count']))
