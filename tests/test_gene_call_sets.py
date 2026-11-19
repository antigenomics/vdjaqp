"""Gene normalization preserves native call ambiguity without inventing ancestry."""
import polars as pl
from vdjaqp.data import to_gene_level,isotype_of


def test_gene_call_sets_retain_alternatives_and_collapse_alleles():
    data=pl.DataFrame({'v_name':['IGHV4-30-4*01,IGHV4-30-8*01','IGHV1-2*01,IGHV1-2*02',None],
                       'j_name':['IGHJ4*01']*3,'c_name':['IGHG1*01']*3})
    result=to_gene_level(data)
    assert result['v_name'].to_list()==['IGHV4-30-4,IGHV4-30-8','IGHV1-2',None]
    assert result['j_name'].to_list()==['IGHJ4']*3


def test_isotype_requires_agreement_across_call_alternatives():
    data=pl.DataFrame({'c_name':['IGHG1*01,IGHG2*01','IGHG1*01,IGHM*01','IGHM*01,IGHM*02',None]})
    assert data.select(isotype_of(data))['isotype'].to_list()==['IgG','unknown','IgM','unknown']
