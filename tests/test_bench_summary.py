import numpy as np

from vdjaqp.bench import summarise


def test_rank_summary_preserves_fractional_median():
    assert summarise(np.array([73, 182]), 33052) == {
        "n_clonotypes": 33052,
        "median_rank": 127.5,
        "worst_rank": 182,
        "top_50": 0,
        "top_500": 2,
    }
    assert summarise(np.array([1, 7, 66]), 33052)["median_rank"] == 7.0


def test_tie_aware_average_precision_auc_and_fixed_shortlist():
    import polars as pl
    import pytest
    from vdjaqp.bench import classification_metrics
    table=pl.DataFrame(dict(clone_id=['a','b','c','d'],score=[2.,2.,1.,0.],confirmed_nats=[20.,21.,19.,0.]))
    m=classification_metrics(table,{'a','c'})
    assert m['average_precision']==pytest.approx(7/12)
    assert m['roc_auc']==pytest.approx(5/8)
    assert (m['true_positive'],m['false_positive'],m['false_negative'],m['true_negative'])==(1,1,1,1)
    assert m['precision_at20']==m['recall_at20']==m['f1_at20']==.5
    assert m['median_rank']==2.
    assert classification_metrics(table.reverse(),{'a','c'})==m
    tied=table.with_columns(pl.lit(0.).alias('score'))
    assert classification_metrics(tied,{'a','c'})['average_precision']==.5
    assert classification_metrics(tied,{'a','c'})['roc_auc']==.5
    with pytest.raises(ValueError):classification_metrics(table,{'absent'})
