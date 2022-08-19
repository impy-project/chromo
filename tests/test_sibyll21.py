from impy.kinematics import EventKinematics
from impy.models import Sibyll21
from impy.constants import TeV
from numpy.testing import assert_allclose, assert_equal
from .util import reference_charge

ekin = EventKinematics(ecm=10 * TeV, p1pdg=2212, p2pdg=2212)
m = Sibyll21(ekin, seed=1)
for event in m(1):
    pass


def test_charge():
    expected = reference_charge(event.id)
    assert_allclose(event.charge, expected)


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


def test_children():
    assert event.children is None


def test_parents():
    assert event.parents is None
