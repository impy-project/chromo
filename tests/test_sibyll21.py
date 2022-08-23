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


def test_to_hepmc3(event):
    hev = event.to_hepmc3()

    assert len(hev.particles) == len(event)
    # no history in Sibyll21
    assert len(hev.vertices) == 0

    for i, p in enumerate(hev.particles):
        assert p.momentum.x == event.px[i]
        assert p.momentum.y == event.py[i]
        assert p.momentum.z == event.pz[i]
        assert p.momentum.e == event.en[i]
        assert p.status == event.status[i]
        assert p.pid == event.pid[i]
        assert p.id == i + 1
