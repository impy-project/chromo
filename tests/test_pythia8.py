from impy.kinematics import CenterOfMass
from impy.models import Pythia8
from impy.constants import GeV, TeV
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from .util import reference_charge
import pytest
from particle import literals as lp


def run_event():
    ekin = CenterOfMass(1 * TeV, 2212, 2212)
    m = Pythia8(ekin, seed=4)
    m.set_stable(lp.pi_0.pdgid, False)
    for event in m(1):
        pass
    return event


# Pythia8 startup time is very slow when all inelastic processes are enabled
_event = run_event()


@pytest.fixture
def event():
    return _event.copy()


def test_charge(event):
    expected = reference_charge(event.pid)
    assert_equal(event.charge, expected)


def test_vertex(event):
    assert np.sum(event.vt != 0) > 0


def test_children(event):
    assert event.children.shape == (len(event), 2)
    # some particles have no children
    assert sum(x[0] == 0 and x[1] == 0 for x in event.children) > 0

    # no particles have single children (no elastic scattering)
    assert sum(x[0] > 0 and x[1] == 0 for x in event.children) == 0

    # some particles have multiple children
    assert sum(x[0] > 0 and x[1] > 0 for x in event.children) > 0


def test_parents(event):
    assert event.parents.shape == (len(event), 2)
    # same particles have no parents
    assert sum(x[0] == 0 and x[1] == 0 for x in event.parents) > 0

    # most particles have a single parent
    assert sum(x[0] > 0 and x[1] == 0 for x in event.parents) > 0

    # some particles have multiple parents
    assert sum(x[0] > 0 and x[1] > 0 for x in event.parents) > 0
