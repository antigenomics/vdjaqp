"""Pure helper checks; never generate a repertoire."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
import csv
import math
from unittest.mock import patch

spec=importlib.util.spec_from_file_location('calibration',Path(__file__).parents[1]/'scripts/calibrate_repertoire.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class Checks(unittest.TestCase):
    def test_inherited_cpu_budget_requires_monitor(self):
        for inherited in ['', '0', '-1', 'wrong', '2.5']:
            with patch.dict(m.os.environ,{'POLARS_MAX_THREADS':inherited}):
                with self.assertRaisesRegex(ValueError,'profile_build'):m.inherited_cpu_budget(16)
        with patch.dict(m.os.environ,{'POLARS_MAX_THREADS':'16'}):
            self.assertEqual(m.inherited_cpu_budget(32),16)
            self.assertEqual(m.inherited_cpu_budget(8),8)
            for requested in [True,0,-1,1.5]:
                with self.assertRaises(ValueError):m.inherited_cpu_budget(requested)

    def test_integer_allocation_distance_and_brackets(self):
        self.assertEqual(m.largest_remainder({'b':1,'a':1},3),{'a':2,'b':1})
        self.assertEqual(m.wasserstein({1:2},{4:3}),dict(numerator=18,denominator=6,distance_AA=3.))
        self.assertEqual(m.wasserstein({1:1,3:1},{2:2})['distance_AA'],1.)
        rows=[]
        for expected in [0.,.1,.4,.05,.025,.013,.0065,.0033]:
            self.assertEqual(m.next_lambda(rows),expected)
            rows.append({'lambda':expected,'W1':{'distance_AA':expected}})
        interior=[{'lambda':x,'W1':{'distance_AA':0 if x==.25 else 1}} for x in [0.,.1,.25,.4]]
        self.assertEqual(m.next_lambda(interior),.18)
        upper=[{'lambda':x,'W1':{'distance_AA':0 if x==.4 else 1}} for x in [0.,.1,.25,.4]]
        self.assertEqual(m.next_lambda(upper),.33)
        for counts,n in [({},1),({'a':0},1),({'a':-1},1),({'a':1},-1)]:
            with self.assertRaises(ValueError):m.largest_remainder(counts,n)
        with self.assertRaises(ValueError):m.wasserstein({}, {1:1})
    def test_precision_limit_and_available_neighbor(self):
        actual=[(0.,4.796741323648103),(.1,2.8844229217110575),(.4,.7296004842615013),(.25,1.4625282485875706),(.33,1.008315173527038),(.37,.8521004842615012),(.39,.7596004842615012)]
        rows=[{'lambda':x,'W1':{'distance_AA':distance}} for x,distance in actual]
        self.assertIsNone(m.next_lambda(rows))
        # Left adjacent interval is rounded to the best point; right is available.
        interior=[{'lambda':x,'W1':{'distance_AA':0 if x==.4 else 1}} for x in [.39,.4,.8]]
        self.assertEqual(m.next_lambda(interior),.6)
        self.assertEqual([m.next_lambda([], (.4,.8,1.6)),m.next_lambda(rows[:1],(.4,.8,1.6)),m.next_lambda(rows[:2],(.4,.8,1.6))],[.4,.8,1.6])
        for starts in [[],[0,.1],[0,.1,.1],[-.1,0,.1],[0,.1,float('inf')],[0,float('nan'),1],[0,True,1],[.4,.1,.8]]:
            with self.assertRaises(ValueError):m.validate_lambda_starts(starts)

    def test_tilt_preserves_pair_measure_and_support(self):
        initial={'np1_lengths':[[n,n+1.] for n in range(41)],'np2_lengths':[[n,41.-n] for n in range(41)]}
        for lam in [.4,.8,1.6]:
            params=m.tilted_parameters(initial,lam)
            for key in initial:
                weights=[w*math.exp(-lam*n) for n,w in initial[key]];total=sum(weights)
                expected={(n,w/total) for (n,_),w in zip(initial[key],weights)}
                self.assertEqual({tuple(p) for p in params[key]},expected)
                self.assertEqual([p[0] for p in sorted(params[key])],list(range(41)))
                self.assertTrue(all(w>0 for n,w in params[key]))
                self.assertEqual(params[key],sorted(params[key],key=lambda p:(p[1],p[0])))
        with self.assertRaises(ValueError):m.tilted_parameters(initial,1000.)

    @unittest.skipUnless(importlib.util.find_spec('GenAIRR'),'pinned generation environment required')
    def test_all_native_initial_plans_compile_without_draws(self):
        import sys
        sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'))
        from GenAIRR.data import HUMAN_IGH_OGRDB as config
        from GenAIRR import Experiment
        native=Experiment.on('human_igh').recombine()._steps[0]
        initial=dict(v_allele_weights={a.name:1. for group in config.v_alleles.values() for a in group},
                     j_allele_weights={a.name:1. for group in config.j_alleles.values() for a in group},
                     np1_lengths=native.np1_lengths,np2_lengths=native.np2_lengths)
        compiled_type=type(Experiment.on('human_igh').recombine().productive_only().compile())
        with patch.object(compiled_type,'run',side_effect=AssertionError('draw prohibited')):
            result=m.preflight_plans(initial,[.4,.8,1.6])
        self.assertEqual([r['lambda'] for r in result],[.4,.8,1.6])
        self.assertTrue(all(r['status']=='compiled' and len(r['parameter_sha256'])==64 for r in result))

    def test_allocation_labels_keep_missing_and_independence(self):
        target={'IGHM':28800,'IGHG1':3882,'IGHG3':223,'IGHA1':64,'IGHE':1,'missing':70}
        a,labels=m.constant_inputs(target,20261111)
        self.assertEqual(len(a),1200);self.assertEqual(len(labels),1200)
        self.assertEqual(sum(x['imputed'] for x in labels),3)
        self.assertEqual(sum(x['observed_C']=='missing' for x in labels),3)
        self.assertEqual((a,labels),m.constant_inputs(target,20261111))
        self.assertNotEqual(list(a.values()),sorted(a.values()))
    def test_explicit_full_constant_inputs(self):
        from collections import Counter
        target={'IGHM':28800,'IGHG1':3882,'IGHG3':223,'IGHA1':64,'IGHE':1,'missing':70}
        default_ids=[f'family{f}_node{n}' for f in range(800) for n in range(1 if f<600 else 3)]
        self.assertEqual(m.constant_inputs(target,20261111),m.constant_inputs(target,20261111,default_ids))
        ids=[f'full_{i}' for i in range(33040)]
        allocation,labels=m.constant_inputs(target,20261112,ids)
        self.assertEqual(set(allocation),set(ids))
        self.assertEqual(Counter(r['observed_C'] for r in labels),Counter(target))
        self.assertEqual(sum(r['imputed'] for r in labels),70)
        self.assertEqual(Counter(r['physical_C'] for r in labels if r['imputed']),Counter(m.largest_remainder({k:v for k,v in target.items() if k!='missing'},70)))
        for invalid in [[],['same','same'],[''],[None]]:
            with self.assertRaises(ValueError):m.constant_inputs(target,20261112,invalid)

    def test_native_event_counts_do_not_sum_rates_or_missing_values(self):
        result=m.native_counter_summary([
            {'native_counters':{'n_mutations':2,'mutation_rate':.01,'n_bad':True}},
            {'native_counters':{'n_mutations':0,'mutation_rate':.02,'n_bad':None}},
            {'native_counters':{'n_mutations':None}},
            {'native_counters':{'n_mutations':1.5}},
            {'native_counters':{}}])
        self.assertNotIn('mutation_rate',result)
        self.assertEqual(result['n_mutations'],dict(sum=2,n_valid=2,n_missing=2,n_invalid=1))
        self.assertEqual(result['n_bad'],dict(sum=None,n_valid=0,n_missing=4,n_invalid=1))

    def test_malformed_targets_fail(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'targets.tsv'
            p.write_text('condition\tpopulation\tquantity\tcategory\tnumerator\tdenominator\tn_observations\nx\tbackground\tV_gene_call\tV\t33039\t33039\t33039\n')
            with self.assertRaises(ValueError):m.targets(p,'x','background')
            with self.assertRaises(ValueError):m.targets(p,'absent','background')

if __name__=='__main__':unittest.main()
