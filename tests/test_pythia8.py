from impy.kinematics import CenterOfMass
from impy.models import Pythia8
from impy.constants import GeV
import numpy as np
from numpy.testing import assert_allclose
from impy.util import AZ2pdg
from .util import reference_charge, run_in_separate_process
import pytest
from particle import literals as lp
from functools import lru_cache


def run_pp_collision():
    evt_kin = CenterOfMass(10 * GeV, "proton", "proton")
    m = Pythia8(evt_kin, seed=4)
    m.set_stable(lp.pi_0.pdgid, True)
    for event in m(1):
        pass

    # some methods only work on original event
    assert np.isnan(event.impact_parameter)
    assert event.n_wounded_A == 0
    assert event.n_wounded_B == 0
    assert event.n_wounded == 0

    return event


@pytest.fixture
@lru_cache(maxsize=1)  # Pythia8 initialization is very slow
def event():
    return run_in_separate_process(run_pp_collision)


def test_charge(event):
    expected = reference_charge(event.pid)
    # skip internal particles unknown to reference_charge
    ma = np.isnan(expected)
    assert np.mean(ma) < 0.1
    event.charge[ma] = np.nan
    assert_allclose(event.charge, expected)


def test_vertex(event):
    assert np.sum(event.vt != 0) > 0


def test_children(event):
    assert event.children.shape == (len(event), 2)
    # some particles have no children
    assert sum(x[0] == 0 and x[1] == 0 for x in event.children) > 0

    # somes particles have single children (elastic scattering in parton shower)
    assert sum(x[0] > 0 and x[1] == 0 for x in event.children) > 0

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


def run_pythia_with_nuclei():
    evt_kin = CenterOfMass(10 * GeV, (4, 2), (12, 6))
    m = Pythia8(evt_kin, seed=1)
    for event in m(1):
        pass
    assert event.impact_parameter > 0
    assert event.n_wounded_A > 0
    assert event.n_wounded_B > 0
    assert event.n_wounded == event.n_wounded_A + event.n_wounded_B


@pytest.mark.skip(reason="it takes forever to generate a nuclear collision")
def test_nuclear_collision():
    run_in_separate_process(run_pythia_with_nuclei)