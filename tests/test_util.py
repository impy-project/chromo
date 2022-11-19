from impy import util
import numpy as np
from numpy.testing import assert_equal


def test_select_parents():
    mask = np.array([False, True, False])
    assert util.select_parents(mask, None) is None

    parents = np.array([[0, 0], [1, 2], [2, 0]])
    par = util.select_parents(mask, parents)
    assert_equal(par, [[0, 0]])

    mask = np.array([False, True, True])
    parents = np.array([[0, 0], [1, 2], [2, 0]])
    par = util.select_parents(mask, parents)
    assert_equal(par, [[0, 0], [1, 0]])


def test_tolerant_string_match():
    assert util.tolerant_string_match("123", "abc1 s2.d34")
    assert not util.tolerant_string_match("123", "321")
    assert not util.tolerant_string_match("123", "3123")
    assert not util.tolerant_string_match("123", "124")
