"""Audit pinned natural constant CDS references; never generate a repertoire.

Raw downloads are inputs. This entry point neither downloads nor modifies them.
"""
import argparse
import hashlib
import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SimpleLocation


MANIFEST = Path(__file__).resolve().parents[1] / 'envs/generation/constant_candidates.json'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def join_constant(variable_nt, variable_aa, constant):
    """Validate native coding variable plus trailing G, without changing any bases."""
    gene = constant['gene']
    require(variable_nt.endswith('G') and len(variable_nt) % 3 == 1,
            f'{gene}: incompatible native variable terminal G/frame')
    joined = variable_nt + constant['nt']
    require(str(Seq(joined).translate()) == variable_aa + constant['aa'] + '*',
            f'{gene}: joined coding translation mismatch')
    return joined


def validated_source_manifest(reference_dir, manifest=MANIFEST):
    """Keep byte pins for sequences; validate the two Ensembl JSON semantics.

    Original JSON bytes are unavailable, so differing bytes are not asserted to
    be key-order-only drift. Every used biological invariant remains pinned.
    """
    data = json.loads(Path(manifest).read_text())
    row = next(r for r in data['candidates'] if r['gene'] == 'IGHM')
    transcript, version = row['nucleotide_accession'].rsplit('.', 1)
    for name, source in data['downloaded_records'].items():
        content = (Path(reference_dir) / name).read_bytes()
        actual = hashlib.sha256(content).hexdigest()
        if actual == source['sha256']:
            continue
        if name == transcript + '.lookup.json':
            meta = json.loads(content)
            expected = dict(id=transcript, version=int(version), assembly_name=row['assembly'],
                            strand=-1, seq_region_name='14', source='havana')
            require(all(meta.get(k) == v for k,v in expected.items()),
                    'IGHM: changed transcript metadata invariants')
            semantic = {k:meta[k] for k in expected}
        elif name == transcript + '.cds-features.json':
            semantic = sorted((f for f in json.loads(content) if f.get('Parent') == transcript),
                              key=lambda f:f['start'], reverse=True)
            require(semantic == row['cds_features'], 'IGHM: changed CDS annotation')
        else:
            raise ValueError(f'{name}: raw SHA256 mismatch')
        source.update(original_sha256=source['sha256'], sha256=actual,
                      validation='pinned biological JSON invariants; full CDS/protein hashes checked downstream',
                      semantic_sha256=hashlib.sha256(json.dumps(semantic,sort_keys=True).encode()).hexdigest())
    return data


def load_validated_constants(reference_dir, manifest=MANIFEST):
    """Extract source records once, verifying raw bytes, annotations and native joins."""
    reference_dir = Path(reference_dir)
    data = validated_source_manifest(reference_dir, manifest)
    for name, source in data['downloaded_records'].items():
        actual = hashlib.sha256((reference_dir / name).read_bytes()).hexdigest()
        require(actual == source['sha256'], f'{name}: raw SHA256 mismatch')
    records = {r.id: r for name in data['downloaded_records'] if name.endswith('.gb')
               for r in SeqIO.parse(reference_dir / name, 'genbank')}
    rows = data['candidates']
    require(len(rows) == 4 and {r['gene'] for r in rows} == {'IGHM', 'IGHG3', 'IGHA1', 'IGHE'},
            'expected exactly four constant candidates')
    result = {}
    for row in rows:
        gene = row['gene']
        if gene == 'IGHM':
            transcript = row['nucleotide_accession'].split('.')[0]
            meta = json.loads((reference_dir / f'{transcript}.lookup.json').read_text())
            require(f"{meta['id']}.{meta['version']}" == row['nucleotide_accession'],
                    'IGHM: transcript version mismatch')
            require(meta['assembly_name'] == row['assembly'], 'IGHM: assembly mismatch')
            features = [f for f in json.loads((reference_dir / f'{transcript}.cds-features.json').read_text())
                        if f.get('Parent') == transcript]
            features.sort(key=lambda f: f['start'], reverse=True)
            require(features == row['cds_features'], 'IGHM: CDS annotation mismatch')
            require(all(f['strand'] == -1 and f['seq_region_name'] == '14' for f in features),
                    'IGHM: expected chromosome 14 reverse strand')
            base = int(row['genomic_region'].split(':')[1].split('..')[0])
            location = CompoundLocation([SimpleLocation(f['start'] - base, f['end'] - base + 1,
                                                        strand=-1) for f in features])
            genomic = SeqIO.read(reference_dir / 'IGHM.GRCh38.genomic.fasta', 'fasta')
            nt = str(location.extract(genomic.seq))
            cds = SeqIO.read(reference_dir / f'{transcript}.cds.fasta', 'fasta')
            require(cds.id == row['nucleotide_accession'] and nt == str(cds.seq),
                    'IGHM: genomic extraction differs from transcript CDS')
            annotated = str(SeqIO.read(reference_dir / (row['uniprot_secreted_isoform'] + '.fasta'), 'fasta').seq)
            offset = 0
        else:
            record = records[row['nucleotide_accession']]
            features = [f for f in record.features if f.type == 'CDS'
                        and f.qualifiers.get('protein_id') == [row['protein_id']]]
            require(len(features) == 1, f'{gene}: expected one matching CDS')
            feature = features[0]
            location = feature.location
            nt = str(feature.extract(record.seq))
            offset = int(feature.qualifiers.get('codon_start', ['1'])[0]) - 1
            annotated = feature.qualifiers['translation'][0]
        require(str(location) == row['location'], f'{gene}: location mismatch')
        require(offset + 1 == row['codon_start'] and len(nt) == row['extracted_nt'],
                f'{gene}: CDS frame/length mismatch')
        require(hashlib.sha256(nt.encode()).hexdigest() == row['cds_sha256'], f'{gene}: CDS hash mismatch')
        require(str(Seq(nt[offset:]).translate()) == annotated + '*', f'{gene}: annotation translation mismatch')
        require(len(annotated) == row['annotated_protein_aa'], f'{gene}: annotation length mismatch')
        if row['native_join'] == 'overlap_one_identical_G':
            require(offset == 0 and nt.startswith('G'), f'{gene}: incompatible overlapping G')
            physical = annotated
        else:
            require(row['native_join'] == 'complete_split_initial_codon' and offset == 2,
                    f'{gene}: unsupported join')
            physical = str(Seq('G' + nt).translate())[:-1]
            require(physical[1:] == annotated, f'{gene}: split codon mismatch')
        require(len(physical) == row['physical_constant_aa'] and
                hashlib.sha256(physical.encode()).hexdigest() == row['physical_constant_sha256'],
                f'{gene}: physical constant hash/length mismatch')
        suffix = nt[1:] if row['native_join'] == 'overlap_one_identical_G' else nt
        require(nt[-3:] in {'TAA', 'TAG', 'TGA'}, f'{gene}: missing authentic terminal stop')
        constant = dict(
            gene=gene, aa=physical, raw_cds_nt=nt, nt=suffix,
            stop_codon=nt[-3:], cds_accession=row['nucleotide_accession'],
            protein_accession=row['protein_id'], cds_sha256=row['cds_sha256'],
            physical_constant_sha256=row['physical_constant_sha256'],
            cds_location=str(location), native_join=row['native_join'])
        # A fixed coding codon plus native trailing G; no generated fixture or base repair.
        join_constant('ATGG', 'M', constant)
        uniprot = str(SeqIO.read(reference_dir / (row['uniprot_secreted_isoform'] + '.fasta'), 'fasta').seq)
        require((physical == uniprot) == row['exact_uniprot_secreted_match'],
                f'{gene}: UniProt comparison mismatch')
        result[gene] = constant
    return result


def validate(reference_dir, manifest=MANIFEST):
    """Return compact verified identities, preserving the audit CLI schema."""
    return [{'gene': r['gene'], 'cds_sha256': r['cds_sha256'],
             'physical_constant_aa': len(r['aa']), 'native_join': r['native_join']}
            for r in load_validated_constants(reference_dir, manifest).values()]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-dir', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(validate(args.reference_dir), indent=2))
