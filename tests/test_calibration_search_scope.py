"""Calibration searches only self-reference texts; its composition prior stays full."""
import numpy as np
import polars as pl
from vdjaqp import errormodel as m
from vdjaqp.alphabet import encode_frame,N_AA
from vdjaqp.match import PeptideMatcher


def test_reference_search_is_exact_with_heldout_db_and_repeated_offsets(monkeypatch):
    clones=pl.DataFrame({'clone_id':['reference','heldout','background'],'seq':['ACDEFGHACDEFGH','QRSTVWYACDEFG','ACDEFGHIKLMNPQRSTVWYYYYY']})
    peptides=pl.DataFrame({'peptide_id':['n','d','h'],'seq':['ACDEFH','ACDEYG','QRSTVW'],'local_confidence':[[.9]*6,[1.]*6,[.9]*6],'source_file':['ref-TRYP.raw','ref-TRYP.raw','hold-TRYP.raw'],'search_type':['denovo','db','denovo']})
    targets=pl.DataFrame({'ms_stem':['ref','hold'],'clone_id':['reference','heldout']})
    monkeypatch.setattr(m,'load_clonotypes',lambda *a:clones)
    monkeypatch.setattr(m,'load_peptides',lambda *a:peptides)
    monkeypatch.setattr(m,'load_targets',lambda *a:targets)
    native=PeptideMatcher.from_clonotypes.__func__;indexed=[];force_full=[True]
    def matcher(cls,frame,k=4):
        indexed.append(frame['clone_id'].to_list())
        return native(cls,clones if force_full[0] else frame,k=k)
    monkeypatch.setattr(PeptideMatcher,'from_clonotypes',classmethod(matcher))
    prior=m.mass_prior;compositions=[]
    def record_prior(composition):
        compositions.append(composition.copy());return prior(composition)
    monkeypatch.setattr(m,'mass_prior',record_prior)
    full=m.ErrorModel.calibrate(exclude={'heldout'})
    force_full[0]=False
    restricted=m.ErrorModel.calibrate(exclude={'heldout'})
    assert indexed==[['reference'],['reference']]
    for field in ['_conf_x','_conf_e','substitution']:
        np.testing.assert_array_equal(getattr(full,field),getattr(restricted,field))
    assert full.n_substitutions==restricted.n_substitutions and full.e_db==restricted.e_db
    code=encode_frame(clones['seq'])[0];expected=np.bincount(code[code>=0],minlength=N_AA).astype(float);expected/=expected.sum()
    for composition in compositions:np.testing.assert_array_equal(composition,expected)
