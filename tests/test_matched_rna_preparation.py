"""Tiny transcript adapter checks; no cohort conversion or reads."""
import copy,hashlib,importlib.util,csv,json,tempfile
from pathlib import Path
import unittest

if __name__!='__main__':
    import pytest
    pytest.importorskip('Bio',reason='pinned generation environment required')

spec=importlib.util.spec_from_file_location('matched_rna',Path(__file__).parents[1]/'scripts/prepare_matched_rna.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class Checks(unittest.TestCase):
    def test_context_model_pin_and_complete_native_coverage(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);source=root/'sources.json';table=root/'contexts.tsv'
            source.write_text(json.dumps(dict(native=[dict(native_allele='V1',locus='IGH')],sources=[dict(source_id='source1',locus='IGH',status='valid',prefix_nt='ATG',leader_nt='ATG')])))
            row=dict(native_allele='V1',locus='IGH',method='STRICT_EXACT',confidence='validated_unique_context',selected_source_id='source1',prefix_nt='ATG',leader_nt='ATG',prefix_sha256=hashlib.sha256(b'ATG').hexdigest())
            with table.open('w') as stream:
                w=csv.DictWriter(stream,fieldnames=list(row),delimiter='\t',lineterminator='\n');w.writeheader();w.writerow(row)
            (root/'audit.json').write_text(json.dumps(dict(model_source_sha256=m.sha(source),native_alleles=1,final_methods={'STRICT_EXACT':1})))
            self.assertEqual(m.load_contexts(table,m.sha(table),source),{'V1':row})
            with self.assertRaisesRegex(ValueError,'hash'):m.load_contexts(table,'wrong',source)
            with table.open('a') as stream:
                w=csv.DictWriter(stream,fieldnames=list(row),delimiter='\t',lineterminator='\n');w.writerow(row)
            with self.assertRaisesRegex(ValueError,'population'):m.load_contexts(table,m.sha(table),source)

    def test_unchanged_mature_substrate_and_truth_offsets(self):
        context=dict(locus='IGH',method='STRICT_EXACT',confidence='validated_unique_context',selected_source_id='public:V1',leader_nt='ATG',prefix_nt='ATG',prefix_sha256=hashlib.sha256(b'ATG').hexdigest())
        row=dict(chain_id='model_H',variable_nt='GCTG',rnaseq_nt='GCTGCTTAA',protein_aa='AA',airr_truth=dict(v_call='V1',v_trim_5=0))
        original=copy.deepcopy(row);seq,provenance=m.convert_chain(row,{'V1':context})
        self.assertEqual(seq,'ATGGCTGCTTAA');self.assertEqual(row,original)
        self.assertEqual((provenance['mature_nt_start'],provenance['leader_nt_start'],provenance['mature_nt_length']),(4,1,9))
        self.assertEqual(provenance['prefix_sha256'],context['prefix_sha256'])
        for field,value in [('variable_nt','ACTG'),('rnaseq_nt','GCTG'),('protein_aa','WRONG')]:
            bad=copy.deepcopy(row);bad[field]=value
            with self.assertRaises(ValueError):m.convert_chain(bad,{'V1':context})
        bad=copy.deepcopy(row);bad['airr_truth']['v_trim_5']=1
        with self.assertRaises(ValueError):m.convert_chain(bad,{'V1':context})

if __name__=='__main__':unittest.main()
