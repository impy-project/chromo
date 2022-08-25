from impy.kinematics import EventKinematics
from impy.models import Pythia6
from impy.constants import GeV, TeV
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from .util import reference_charge, run_in_separate_process
import pytest
import pickle
from particle import literals as lp


def run_event():
    ekin = EventKinematics(ecm=1 * TeV, p1pdg=2212, p2pdg=2212)
    m = Pythia6(ekin, seed=4)
    m.set_stable(lp.pi_0.pdgid, False)  # needed to get nonzero vertices
    for event in m(100):
        if len(event) > 10:  # to skip elastic events
            break
    return event  # MCEvent is pickeable, but restored as EventData


@pytest.fixture
def event():
    return run_in_separate_process(run_event)


def test_charge(event):
    expected = reference_charge(event.pid)
    assert_allclose(event.charge, expected)


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


def run_is_view():
    ekin = EventKinematics(ecm=10 * GeV, p1pdg=2212, p2pdg=2212)

    m = Pythia6(ekin, seed=1)
    for event in m(1):
        pass

    return (
        event.px.flags["OWNDATA"],
        event[:5].px.flags["OWNDATA"],
        event.copy().px.flags["OWNDATA"],
    )


def test_is_view():
    (
        event_owndata,
        sliced_owndata,
        copy_owndata,
    ) = run_in_separate_process(run_is_view)
    assert event_owndata is False
    assert sliced_owndata is False
    assert copy_owndata is True


def test_final_state(event):
    ev1 = event.final_state()
    ev2 = event[event.status == 1]
    assert_equal(ev1, ev2)


def test_final_state_charged(event):
    ev1 = event.final_state_charged()
    ev2 = event[(event.status == 1) & (event.charge != 0)]
    ev3 = event[event.status == 1]
    ev3 = ev3[ev3.charge != 0]
    assert_equal(ev1, ev2)
    assert_equal(ev1, ev3)


def run_pickle():
    ekin = EventKinematics(ecm=10 * GeV, p1pdg=2212, p2pdg=2212)

    m = Pythia6(ekin, seed=1)
    for event in m(1):
        pass

    s = pickle.dumps(event)
    event2 = pickle.loads(s)

    assert event == event2


def test_pickle(event):
    run_in_separate_process(run_pickle)
