"""Reference annotation tests: metadata only, never native cohort generation."""
import csv
from pathlib import Path
import sys
import tempfile
import unittest
from importlib.util import find_spec

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from genairr_gene_map import (exact_v_map, canonical_v_fields,
                             load_reference_gene_mapping, write_reference_gene_mapping)


class GeneMapTests(unittest.TestCase):
    def test_exact_join_and_scalars(self):
        legacy = [('L1', 'group', 'aaa'), ('L2', 'group', 'CCC'),
                  ('L3', 'group', 'GGG'), ('L4', 'group', 'TTT')]
        canonical = [('A*1', 'A', 'AAA'), ('B*1', 'B', 'CCC'), ('B*2', 'B', 'CCC'),
                     ('C*1', 'C', 'GGG'), ('D*1', 'D', 'GGG')]
        rows = exact_v_map(legacy, canonical)
        self.assertEqual([r['canonical_v_match_status'] for r in rows.values()],
                         ['unique_gene', 'unique_gene', 'ambiguous_genes', 'no_exact_match'])
        self.assertEqual(rows['L2']['canonical_alleles'], ['B*1', 'B*2'])
        self.assertEqual(rows['L2']['canonical_v_gene'], 'B')
        self.assertEqual(rows['L3']['canonical_genes'], ['C', 'D'])
        self.assertIsNone(rows['L3']['canonical_v_gene'])
        self.assertEqual(rows['L4']['canonical_v_gene_count'], 0)
        mapping = {'alleles': rows}
        self.assertEqual(canonical_v_fields('IGH', 'L2', mapping)['canonical_v_gene_count'], 1)
        self.assertEqual(canonical_v_fields('IGK', 'anything', mapping),
                         dict(canonical_v_gene=None, canonical_v_gene_count=None,
                              canonical_v_match_status='not_evaluated_locus'))
        with self.assertRaises(KeyError):
            canonical_v_fields('IGH', 'unknown', mapping)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'map.tsv'
            write_reference_gene_mapping(path, mapping)
            with path.open() as stream:
                flat = list(csv.DictReader(stream, delimiter='\t'))
            self.assertEqual(len(flat), 6)
            self.assertEqual(flat[-1]['canonical_allele'], 'null')
            self.assertEqual(flat[-1]['canonical_gene'], 'null')

    @unittest.skipUnless(find_spec('GenAIRR'), 'requires pinned generation environment')
    def test_installed_metadata(self):
        mapping = load_reference_gene_mapping()
        self.assertEqual(mapping['legacy_allele_count'], 198)
        self.assertEqual(mapping['canonical_allele_count'], 190)
        statuses = [r['canonical_v_match_status'] for r in mapping['alleles'].values()]
        self.assertEqual([statuses.count(x) for x in
                          ('unique_gene', 'ambiguous_genes', 'no_exact_match')], [163, 7, 28])
        identity = Path(__file__).resolve().parents[1]/'envs/generation/identity.json'
        import json
        bad = json.loads(identity.read_text())
        bad['canonical_gene_reference_sha256']['human_igh_ogrdb_10_observed_v1.zip'] = '0'*64
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'identity.json'
            path.write_text(json.dumps(bad))
            with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                load_reference_gene_mapping(path)


if __name__ == '__main__':
    unittest.main()
