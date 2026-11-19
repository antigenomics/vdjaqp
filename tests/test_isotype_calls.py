"""Producer-level class calls and missing C evidence retain their semantics."""
import polars as pl
from vdjaqp.data import isotype_of


def test_generic_and_subclass_igg_calls():
    frame=pl.DataFrame({'c_name':['IGHG','IGHG1*01','IGHG4','IGHM','IGHA1',None,'unknown','IGHC']})
    assert frame.select(isotype_of(frame))['isotype'].to_list()==['IgG','IgG','IgG','IgM','other','unknown','unknown','unknown']


def test_missing_isotype_is_neutral_and_does_not_change_observed_odds():
    import numpy as np
    from vdjaqp.aggregate import isotype_log_prior
    known=np.array(['IgG','IgM','other'])
    original=np.where(known=='IgG',np.log((0.9/1)/((1-0.9)/2)),0.0)
    np.testing.assert_array_equal(isotype_log_prior(known),original)
    missing=isotype_log_prior(np.array(['IgG','IgM','other','unknown','unknown']))
    np.testing.assert_array_equal(missing[-2:],[0.,0.])
    assert np.isclose(missing[0]-missing[1],original[0]-original[1])
    np.testing.assert_array_equal(isotype_log_prior(np.array(['IgG','unknown'])),[0.,0.])
