from impy.kinematics import EventKinematics
from impy.models import Sibyll21
from impy.constants import TeV
from numpy.testing import assert_allclose, assert_equal
from .util import reference_charge, run_in_separate_process
import pytest


def event_run():
    ekin = EventKinematics(ecm=10 * TeV, p1pdg=2212, p2pdg=2212)
    m = Sibyll21(ekin, seed=1)
    for event in m(1):
        pass
    return event.copy()  # copy is pickleable


@pytest.fixture
def event():
    return run_in_separate_process(event_run)


def test_charge(event):
    expected = reference_charge(event.pid)
    assert_allclose(event.charge, expected)


def test_children(event):
    assert event.children is None


def test_parents(event):
    assert event.parents is None


def test_vertex(event):
    # no vertex info in Sibyll21
    assert_equal(event.vx, 0)
    assert_equal(event.vy, 0)
    assert_equal(event.vz, 0)
    assert_equal(event.vt, 0)
