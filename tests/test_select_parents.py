from impy.util import select_parents
import numpy as np
from numpy.testing import assert_equal


def test_select_parents():
    mask = np.array([False, True, False])
    assert select_parents(mask, None) is None

    parents = np.array([[0, 0], [1, 2], [2, 0]])
    par = select_parents(mask, parents)
    assert_equal(par, [[0, 0]])

    mask = np.array([False, True, True])
    parents = np.array([[0, 0], [1, 2], [2, 0]])
    par = select_parents(mask, parents)
    assert_equal(par, [[0, 0], [1, 0]])
