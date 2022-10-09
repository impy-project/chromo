from impy.kinematics import CenterOfMass
from impy.models import Pythia8
from impy.constants import GeV, TeV
import numpy as np
from numpy.testing import assert_allclose
from .util import reference_charge, run_in_separate_process
import pytest
from particle import literals as lp
from functools import cache


def run_event():
    ekin = CenterOfMass(1 * TeV, 2212, 2212)
    m = Pythia8(ekin, seed=4)
    m.set_stable(lp.pi_0.pdgid, False)
    for event in m(1):
        pass
    return event


@pytest.fixture
@cache  # Pythia8 initialization is very slow
def event():
    return run_in_separate_process(run_event)


def test_charge(event):
    expected = reference_charge(event.pid)
    # skip internal particles unknown to reference_charge
    ma = np.isnan(expected)
    assert np.mean(ma) < 0.1
    event.charge[ma] = np.nan
    assert_allclose(event.charge, expected)


def test_vertex(event):
    assert np.sum(event.vt != 0) > 0


def test_parents(event):
    assert event.parents.shape == (len(event), 2)
    # same particles have no parents
    assert sum(x[0] == 0 and x[1] == 0 for x in event.parents) > 0

    # most particles have a single parent
    assert sum(x[0] > 0 and x[1] == 0 for x in event.parents) > 0

    # some particles have multiple parents
    assert sum(x[0] > 0 and x[1] > 0 for x in event.parents) > 0
