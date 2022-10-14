from impy.kinematics import CenterOfMass
from impy.models import EposLHC
from impy.constants import GeV
import numpy as np
from numpy.testing import assert_allclose
from .util import reference_charge, run_in_separate_process
import pytest
from particle import literals as lp
from functools import lru_cache


def run_pp_collision():
    evt_kin = CenterOfMass(10 * GeV, "proton", "proton")
    m = EposLHC(evt_kin, seed=4)
    m.set_stable(lp.pi_0.pdgid, True)
    for event in m(1):
        pass

    # some methods only work on original event
    assert event.impact_parameter > 0
    assert event.n_wounded_A == 1
    assert event.n_wounded_B == 1
    assert event.n_wounded == 2

    return event


@pytest.fixture
@lru_cache(maxsize=1)  # EposLHC initialization is very slow
def event():
    return run_in_separate_process(run_pp_collision)


def test_charge(event):
    expected = reference_charge(event.pid)
    # skip internal particles unknown to reference_charge
    ma = np.isnan(expected)
    assert np.mean(ma) < 0.3
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


def run_ion_collision():
    evt_kin = CenterOfMass(10 * GeV, (4, 2), (12, 6))
    m = EposLHC(evt_kin, seed=1)
    for event in m(1):
        pass
    assert event.impact_parameter > 0
    assert event.n_wounded_A > 0
    assert event.n_wounded_B > 0
    assert event.n_wounded == event.n_wounded_A + event.n_wounded_B


@pytest.mark.skip(reason="Bug: EposLHC initialization with nuclei fails")
def test_ion_collision():
    run_in_separate_process(run_ion_collision)
