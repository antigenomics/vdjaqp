"""The sampling model preserves fixed depth and never creates unseen identities."""
import importlib.util
from pathlib import Path
import numpy as np

spec = importlib.util.spec_from_file_location('sampling', Path(__file__).parents[1]/'scripts/simulate_repertoire_sampling.py')
sampling = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sampling)


def test_fixed_depth_and_zero_probability_identity():
    values = sampling.draw_library((np.array([7,3,0]),10,100,16))
    assert np.all(values.sum(axis=1)==10)
    assert np.all(values[:,2]==0)
    assert np.unique(values[:,0]).size>1
    assert np.array_equal(values,sampling.draw_library((np.array([7,3,0]),10,100,16)))


def test_unresolved_remainder_is_not_a_junction():
    values = sampling.draw_library((np.array([0,0]),10,100,16))
    assert not values.any()
