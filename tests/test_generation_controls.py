"""No-draw checks: run with .tool-envs/generation/bin/python tests/test_generation_controls.py."""
import copy
from dataclasses import asdict
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

if __name__ != '__main__':
    import pytest
    pytest.importorskip('GenAIRR', reason='pinned generation environment required')

from GenAIRR import Experiment
from GenAIRR.data import HUMAN_IGH_OGRDB as config
repo = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('population_controls', repo/'scripts/generate_population.py')
pop = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pop)


class Controls(unittest.TestCase):
    def test_exact_layout_and_constant_allocation(self):
        layout = pop.mixture_layout(33040, memory_records=16518)
        self.assertEqual(layout['naive_antibodies'], 16522)
        self.assertEqual(layout['families'], 22028)
        self.assertEqual(layout['memory_families'], 5506)
        for kwargs in [dict(memory_records=1), dict(memory_records=True),
                       dict(memory_records=33042), dict(memory_records=12, memory_fraction=.5)]:
            with self.assertRaises(ValueError):pop.mixture_layout(33040, **kwargs)
        small = pop.mixture_layout(24, memory_records=12)
        ids = [f'family{f}_node{n}' for f in range(small['families'])
               for n in range(1 if f < small['naive_antibodies'] else 3)]
        allocation = {name:pop.HEAVY_CONSTANT_GENES[i % 5] for i,name in enumerate(ids)}
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'allocation.json'
            model = dict(schema_version=1, allocation=allocation)
            path.write_text(json.dumps(model))
            loaded, receipt = pop.load_heavy_allocation(path, small)
            self.assertEqual(loaded, allocation)
            self.assertEqual(receipt['sha256'], hashlib.sha256(path.read_bytes()).hexdigest())
            for bad in [dict(schema_version=True, allocation=allocation),
                        dict(schema_version=1, allocation={}),
                        dict(schema_version=1, allocation=dict(allocation, extra='IGHM')),
                        dict(schema_version=1, allocation=dict(allocation, **{ids[0]:'INVALID'}))]:
                path.write_text(json.dumps(bad))
                with self.assertRaises(ValueError):pop.load_heavy_allocation(path, small)
            path.write_text(json.dumps(model).replace('"schema_version": 1', '"schema_version": 1, "schema_version": 1'))
            with self.assertRaisesRegex(ValueError,'duplicate'):pop.load_heavy_allocation(path,small)

    def test_native_model_validation_and_compile_without_draws(self):
        identity = json.loads((repo/'envs/generation/identity.json').read_text())
        baseline_exp = Experiment.on('human_igh').recombine()
        baseline_step = asdict(baseline_exp._steps[0])
        model = dict(schema_version=1, genairr_version=identity['versions']['GenAIRR'],
                     reference_sha256=identity['genairr_reference_sha256']['HUMAN_IGH_OGRDB.pkl'],
                     parameters=dict(v_allele_weights={a.name: 1. for group in config.v_alleles.values() for a in group},
                                     j_allele_weights={a.name: 1. for group in config.j_alleles.values() for a in group},
                                     np1_lengths=[list(pair) for pair in baseline_step['np1_lengths']],
                                     np2_lengths=[list(pair) for pair in baseline_step['np2_lengths']]))
        self.assertEqual(pop.load_recombination_parameters(None), ({}, None))
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'model.json'
            path.write_text(json.dumps(model))
            kwargs, provenance = pop.load_recombination_parameters(path)
            self.assertEqual(kwargs, model['parameters'])
            self.assertEqual(provenance['sha256'], hashlib.sha256(path.read_bytes()).hexdigest())
            # Pinned GenAIRR exposes no public numeric-plan serializer. Inspect its
            # pure-data recombination IR; None usage lowers to native uniform.
            explicit = asdict(Experiment.on('human_igh').recombine(**kwargs)._steps[0])
            expected = dict(baseline_step)
            expected.update(weights_v=tuple([1.]*len(kwargs['v_allele_weights'])),
                            weights_j=tuple([1.]*len(kwargs['j_allele_weights'])),
                            np_lengths_overridden=True)
            self.assertEqual(explicit, expected)
            changed = copy.deepcopy(kwargs)
            changed['v_allele_weights'][next(iter(changed['v_allele_weights']))] = 2.
            changed['np1_lengths'][0][1] *= 2.
            altered = asdict(Experiment.on('human_igh').recombine(**changed)._steps[0])
            self.assertEqual({k for k in explicit if explicit[k] != altered[k]},
                             {'weights_v', 'np1_lengths'})
            self.assertEqual(sum(a != b for a, b in zip(explicit['weights_v'], altered['weights_v'])), 1)
            self.assertEqual(sum(a != b for a, b in zip(explicit['np1_lengths'], altered['np1_lengths'])), 1)
            baseline = baseline_exp.productive_only().compile()
            native = Experiment.on('human_igh').recombine(**kwargs).productive_only().compile()
            self.assertEqual(native.pass_plan, baseline.pass_plan)
            self.assertEqual(native.active_contracts, baseline.active_contracts)
            self.assertEqual(native.simulator.declared_choices(), baseline.simulator.declared_choices())
            duplicate = json.dumps(model).replace('"schema_version": 1', '"schema_version": 1, "schema_version": 1')
            path.write_text(duplicate)
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                pop.load_recombination_parameters(path)
            allele = next(iter(model['parameters']['v_allele_weights']))
            path.write_text(json.dumps(model).replace('"'+allele+'": 1.0', '"'+allele+'": 1.0, "'+allele+'": 2.0'))
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                pop.load_recombination_parameters(path)
            invalid = []
            for field, value in [('schema_version', True), ('genairr_version', 'wrong'), ('reference_sha256', 'wrong')]:
                bad = copy.deepcopy(model); bad[field] = value; invalid.append(bad)
            for key in ['v_allele_weights', 'j_allele_weights']:
                bad = copy.deepcopy(model); bad['parameters'][key].pop(next(iter(bad['parameters'][key]))); invalid.append(bad)
                for weight in [True, 0, -1, float('nan'), float('inf')]:
                    bad = copy.deepcopy(model); bad['parameters'][key][next(iter(bad['parameters'][key]))] = weight; invalid.append(bad)
            for pairs in [[[n, 1.] for n in range(40)], [[0, 1.]]*41, [[n, True] for n in range(41)]]:
                bad = copy.deepcopy(model); bad['parameters']['np1_lengths'] = pairs; invalid.append(bad)
            bad = copy.deepcopy(model); bad['parameters']['d_allele_weights'] = {}; invalid.append(bad)
            for bad in invalid:
                path.write_text(json.dumps(bad))
                with self.assertRaises(ValueError):pop.load_recombination_parameters(path)


if __name__ == '__main__':
    unittest.main()
