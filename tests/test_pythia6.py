from impy.kinematics import EventKinematics
from impy.models import Pythia6
from impy.constants import TeV
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from .util import reference_charge
import pytest
import pickle

ekin = EventKinematics(ecm=10 * TeV, p1pdg=2212, p2pdg=2212)
m = Pythia6(ekin, seed=1)
for event in m(1):
    pass


def test_charge():
    expected = reference_charge(event.pid)
    assert_allclose(event.charge, expected)


@pytest.mark.xfail(reason="no vertex info in pythia6")
def test_vertex():
    assert np.sum(event.vt != 0) > 0


def test_children():
    assert event.children.shape == (len(event), 2)
    # some particles have no children
    assert sum(x[0] == 0 and x[1] == 0 for x in event.children) > 0

    # no particles have single children (no elastic scattering)
    assert sum(x[0] > 0 and x[1] == 0 for x in event.children) == 0

    # some particles have multiple children
    assert sum(x[0] > 0 and x[1] > 0 for x in event.children) > 0


def test_parents():
    assert event.parents.shape == (len(event), 2)
    # same particles have no parents
    assert sum(x[0] == 0 and x[1] == 0 for x in event.parents) > 0

    # most particles have a single parent
    assert sum(x[0] > 0 and x[1] == 0 for x in event.parents) > 0

    # some particles have multiple parents
    assert sum(x[0] > 0 and x[1] > 0 for x in event.parents) > 0


def test_event_is_readonly_view():
    assert event.px.flags["OWNDATA"] is False
    assert event.px.flags["WRITEABLE"] is False


def test_final_state():
    ev1 = event.final_state()
    ev2 = event[event.status == 1]
    assert_equal(ev1, ev2)


def test_final_state_charged():
    ev1 = event.final_state_charged()
    ev2 = event[(event.status == 1) & (event.charge != 0)]
    ev3 = event[event.status == 1]
    ev3 = ev3[ev3.charge != 0]
    assert_equal(ev1, ev2)
    assert_equal(ev1, ev3)


def test_pickle():
    # cannot pickle MCEvent...
    with pytest.raises(TypeError):
        pickle.dumps(event)

    # but can pickle EventData
    ev1 = event.copy()

    s = pickle.dumps(ev1)
    ev2 = pickle.loads(s)

    assert ev1 == ev2
