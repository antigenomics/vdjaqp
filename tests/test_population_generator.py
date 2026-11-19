"""Run with the pinned generation interpreter: python tests/test_population_generator.py."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import json
import shutil
import os
from collections import Counter

if __name__ != '__main__':
    import pytest
    pytest.importorskip('Bio', reason='pinned generation environment required')
    pytest.importorskip('GenAIRR', reason='pinned generation environment required')

from Bio.Seq import translate

spec = importlib.util.spec_from_file_location('population', Path(__file__).parents[1]/'scripts/generate_population.py')
pop = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pop)

class PopulationChecks(unittest.TestCase):
    def test_invalid_never_writes(self):
        for kwargs in [{'families': 3}, {'families': 0}, {'pairs_per_chain': 0}, {'mutation_rate': float('nan')}, {'mutation_rate': -1},
                       {'heavy_constant_gene': 'unknown'},
                       {'heavy_constant_gene': 'IGHM'},
                       {'additional_reference_dir': 'unused'}]:
            with tempfile.TemporaryDirectory() as temp:
                output = Path(temp)/'absent'
                with self.assertRaises(ValueError):
                    pop.generate(output, temp, **dict({'families': 2}, **kwargs))
                self.assertFalse(output.exists())

    def test_reference_hash_failure_never_writes(self):
        with tempfile.TemporaryDirectory() as temp:
            raw = Path(temp)/'raw'
            raw.mkdir()
            for name in ('P01857-1.fasta','P01834.fasta','J00228.embl','J00241.embl','P0DOY2.fasta','J00253.embl'):
                (raw/name).write_text('wrong raw bytes')
            output = Path(temp)/'absent'
            with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                pop.generate(output, raw, families=2)
            self.assertFalse(output.exists())

    def test_constant_selection_without_native_draws(self):
        repo = Path(__file__).resolve().parents[1]
        raw = repo/'tmp/resource-120-complete/references'
        additional = repo/'tmp/constant-reference-audit'
        if not (raw/'J00228.embl').exists() or not additional.exists():
            self.skipTest('pinned raw reference inputs required')
        with patch.object(pop.Experiment, 'on', side_effect=AssertionError('native draw')):
            default, _, extra = pop.load_population_constants(raw)
            self.assertEqual(default, pop.load_constants(raw, all_loci=True))
            self.assertIsNone(extra)
            with tempfile.TemporaryDirectory() as temp:
                target = Path(temp)/'constants.fasta'
                pop.write_constant_references(target, default)
                old_records = []
                for name in ('P01857-1.fasta', 'P01834.fasta', 'P0DOY2.fasta'):
                    lines = (raw/name).read_text().splitlines()
                    accession = lines[0].split('|')[1]
                    old_records.append(f'>constant_{accession}\n'+''.join(lines[1:])+'\n')
                self.assertEqual(target.read_text(), ''.join(old_records))
            validated = pop.load_validated_constants(additional)
            for gene in ('IGHM', 'IGHG3', 'IGHA1', 'IGHE'):
                selected, _, sources = pop.load_population_constants(raw, gene, additional)
                self.assertEqual(selected['IGH'], validated[gene])
                with tempfile.TemporaryDirectory() as temp:
                    target = Path(temp)/'constants.fasta'
                    pop.write_constant_references(target, selected)
                    self.assertEqual(target.read_text().splitlines()[:2],
                                     [f">constant_{validated[gene]['protein_accession']}", validated[gene]['aa']])
                    self.assertEqual(len(target.read_text().splitlines()), 6)
                self.assertEqual(selected['IGK'], default['IGK'])
                self.assertEqual(selected['IGL'], default['IGL'])
                record = {'sequence': 'ATGG', 'sequence_aa': 'M', 'productive': True}
                self.assertIsNone(pop.compatibility_reason(record, selected['IGH']))
                self.assertIsNotNone(pop.compatibility_reason(dict(record, sequence='ATGA'), selected['IGH']))
                self.assertEqual(sources, json.loads(pop.CONSTANT_MANIFEST.read_text()))
            with tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                copied = root/'additional'
                copied.mkdir()
                for name in sources['downloaded_records']:
                    shutil.copyfile(additional/name, copied/name)
                # An unused candidate's raw source must still be checked.
                (copied/'P01854-1.fasta').write_text('tampered')
                with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
                    pop.generate(root/'absent', raw, families=2,
                                 heavy_constant_gene='IGHM', additional_reference_dir=copied)
                self.assertFalse((root/'absent').exists())

    def test_mixture_layout_counts_and_balanced_abundance(self):
        for fraction, naive, trees in ((.25, 900, 100), (.5, 600, 200), (.75, 300, 300)):
            layout = pop.mixture_layout(1200, fraction, 'head80')
            self.assertEqual(layout['naive_antibodies'], naive)
            self.assertEqual(layout['memory_families'], trees)
            self.assertEqual(naive + 3*trees, 1200)
            for start, count, size in ((0, naive, 1), (naive, trees, 3)):
                for light in (0, 1):
                    families = range(start+light, start+count, 2)
                    head = [i for i in families if i in layout['head_families']]
                    self.assertEqual(len(head)*10, len(families))
                    total = sum((8 if i in head else 2/9)*size for i in families)
                    self.assertAlmostEqual(total, len(families)*size)
                    self.assertAlmostEqual(len(head)*8*size/total, .8)
        self.assertEqual(pop.mixture_layout(12, .5)['memory_families'], 2)
        for n, f, abundance in ((12, .5, 'head80'), (10, .25, 'equal'),
                                (0, .5, 'equal'), (12, .4, 'equal')):
            with self.assertRaises(ValueError):
                pop.mixture_layout(n, f, abundance)

    def test_invalid_mixture_never_writes(self):
        cases = [dict(total_antibodies=12, memory_fraction=.5),
                 dict(total_antibodies=12, memory_fraction=.5, families=2, additional_reference_dir='raw'),
                 dict(total_antibodies=12, memory_fraction=.5, protein_concentration='head80', additional_reference_dir='raw'),
                 dict(total_antibodies=12, memory_fraction=.5, heavy_constant_gene='IGHE', additional_reference_dir='raw')]
        for kwargs in cases:
            with tempfile.TemporaryDirectory() as temp:
                with self.assertRaises(ValueError):
                    pop.generate(Path(temp)/'absent', temp, **kwargs)
                self.assertFalse((Path(temp)/'absent').exists())

    @unittest.skipUnless(os.environ.get('VDJAQP_NATIVE_MIXTURE_TEST') == '1',
                         'explicit native engineering fixture only')
    def test_tiny_native_mixture(self):
        repo = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)/'population'
            rows = pop.generate(output, repo/'tmp/resource-120-complete/references',
                total_antibodies=12, memory_fraction=.5,
                additional_reference_dir=repo/'tmp/constant-reference-audit')
            truth = json.loads((output/'truth.json').read_text())
            self.assertEqual(len(rows), 24)
            self.assertEqual(Counter(r['component'] for r in rows), {'naive': 12, 'memory': 12})
            self.assertEqual(len({r['antibody_id'] for r in rows}), 12)
            self.assertEqual(len(truth['latent_ancestors']), 4)
            all_rows = {r['chain_id']: r for r in rows + truth['latent_ancestors']}
            for row in rows:
                self.assertEqual(str(translate(row['rnaseq_nt'])), row['protein_aa']+'*')
                self.assertEqual(row['rna_pairs']*row['secretion_factor'], row['protein_intensity'])
                if row['component'] == 'naive':
                    self.assertIsNone(row['parent_id'])
                    self.assertEqual(row['root_distance']['substitutions'], 0)
                else:
                    ancestor = all_rows[row['ancestor_id']]
                    parent = all_rows[row['parent_id']]
                    self.assertEqual(row['root_distance'], pop.distance(row['variable_nt'], ancestor['variable_nt']))
                    self.assertEqual(row['branch_distance'], pop.distance(row['variable_nt'], parent['variable_nt']))
                if row['locus'] == 'IGH':
                    self.assertEqual(row['constant_nt_gene'], 'IGHM' if row['component'] == 'naive' else 'IGHG1')
            self.assertEqual(len((output/'read_counts.tsv').read_text().splitlines()), 24)
            fasta = (output/'proteins.fasta').read_text()
            self.assertEqual(fasta.count('>'), 24)
            self.assertNotIn('_ancestor_', fasta)
            self.assertEqual((output/'constant_reference_proteins.fasta').read_text().count('>'), 4)
            seeds = [a['seed'] for a in truth['attempts']]
            self.assertEqual(len(seeds), len(set(seeds)))

    def test_boundary_condition(self):
        good = {'sequence': 'GCTG', 'sequence_aa': 'A', 'productive': True}
        constant = {'nt': 'CTTAA', 'aa': 'A'}
        self.assertIsNone(pop.compatibility_reason(good, constant))
        self.assertEqual(pop.compatibility_reason(dict(good, sequence='GCTA'), constant), 'constant_translation_mismatch')

    def test_distance(self):
        self.assertEqual(pop.distance('ACGT', 'AGGT'),
                         {'substitutions': 1, 'denominator_nt': 4, 'fraction': .25})
        self.assertEqual(pop.distance('ACGT', 'ACGT')['substitutions'], 0)
        with self.assertRaises(ValueError):
            pop.distance('ACGT', 'ACG')

    def test_seed_ranges(self):
        seeds = [pop.draw_seed(20261110, stage, locus, family, node, attempt)
                 for stage in (1, 2) for locus in range(3) for family in range(4)
                 for node in range(3) for attempt in range(100)]
        self.assertEqual(len(seeds), len(set(seeds)))
        limits = [10**14, 3, 3, 10**10, 3, 100]
        for index, upper in enumerate(limits):
            for invalid in (-1, upper):
                coords = [20261110, 1, 0, 0, 0, 0]
                coords[index] = invalid
                with self.assertRaises(ValueError):
                    pop.draw_seed(*coords)
        self.assertLess(pop.draw_seed(10**14-1, 2, 2, 10**10-1, 2, 99), 2**64)
        root_max = (10**14-1) + 2*10**14 + (10**10-1)
        self.assertLess(root_max, pop.draw_seed(10**14-1, 1, 0, 0, 0, 0))

if __name__ == '__main__':
    unittest.main()
