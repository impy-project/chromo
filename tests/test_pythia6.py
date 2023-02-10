# These tests check Pythia6 but also the MCEvent in general.
# It is not necessary to duplicate all tests for every model.
from chromo.kinematics import CenterOfMass
from chromo.models import Pythia6
from chromo.constants import GeV, TeV
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from .util import reference_charge
import pytest
import pickle
from functools import lru_cache


def test_name():
    assert Pythia6.name == "Pythia"
    assert Pythia6.version == "6.428"
    assert Pythia6.label == "Pythia-6.428"


def run_collision(p1, p2):
    kin = CenterOfMass(100 * GeV, p1, p2)
    m = Pythia6(seed=4)
    for event in m(kin, 1):
        pass
    return event  # MCEvent is restored as EventData


@pytest.fixture
@lru_cache(maxsize=1)
def event():
    return run_collision("p", "p").copy()


def test_cross_section():
    kin = CenterOfMass(10 * GeV, "p", "p")
    m = Pythia6(seed=1)
    c = m.cross_section(kin)

    assert_allclose(c.total, 38.4, atol=0.1)
    assert_allclose(c.inelastic, 31.4, atol=0.1)
    assert_allclose(c.elastic, 7.0, atol=0.1)
    assert_allclose(c.diffractive_xb, 2.6, atol=0.1)
    assert_allclose(c.diffractive_ax, 2.6, atol=0.1)
    assert_allclose(c.diffractive_xx, 0.9, atol=0.1)
    assert c.diffractive_axb == 0
    assert_allclose(
        c.non_diffractive,
        c.inelastic - c.diffractive_xb - c.diffractive_ax - c.diffractive_xx,
    )


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


def test_is_view():
    kin = CenterOfMass(10 * GeV, 2212, 2212)

    m = Pythia6(seed=1)
    for event in m(kin, 1):
        pass

    event.px.flags["OWNDATA"] is False
    event[:5].px.flags["OWNDATA"] is False
    event.copy().px.flags["OWNDATA"] is True


def test_final_state(event):
    ev1 = event.final_state()
    ev2 = event[event.status == 1]
    ev2.parents = None  # event.final_state() drops parents
    assert_equal(ev1, ev2)


def test_final_state_charged(event):
    ev1 = event.final_state_charged()
    ev2 = event[(event.status == 1) & (event.charge != 0)]
    ev3 = event[event.status == 1]
    ev3 = ev3[ev3.charge != 0]
    ev2.parents = None  # event.final_state() drops parents
    ev3.parents = None  # event.final_state() drops parents
    assert_equal(ev1, ev2)
    assert_equal(ev1, ev3)


def test_pickle():
    kin = CenterOfMass(10 * GeV, 2212, 2212)

    m = Pythia6(seed=1)
    for event in m(kin, 1):
        pass

    s = pickle.dumps(event)

    event1 = event.copy()

    # draw another event so that MCEvent
    # get overridden
    for _ in m(kin, 1):
        pass

    event2 = pickle.loads(s)

    assert event1 == event2


def test_event_getstate_setstate():
    from impy.models.pythia6 import PYTHIA6Event
    from impy.common import EventData

    kin = CenterOfMass(10 * GeV, 2212, 2212)

    m = Pythia6(seed=1)
    for event in m(kin, 1):
        pass

    assert type(event) is PYTHIA6Event

    # this must create a copy
    state = event.__getstate__()

    event1 = event.copy()
    assert type(event1) is not PYTHIA6Event

    # draw another event so that MCEvent gets overridden
    for _ in m(kin, 1):
        pass

    event2 = EventData(**state)

    assert event1 == event2


def test_event_copy():
    from impy.models.pythia6 import PYTHIA6Event
    from impy.common import EventData

    m = Pythia6(seed=4)
    m.set_stable("pi_0", False)  # needed to get nonzero vertices
    kin = CenterOfMass(1 * TeV, 2212, 2212)
    for event in m(kin, 1):
        pass

    event2 = event.copy()
    assert type(event) is PYTHIA6Event
    assert type(event2) is not PYTHIA6Event
    assert type(event2) is EventData
    assert event2 is not event
    assert event2 == event

    event3 = event2.copy()
    assert type(event3) is EventData
    assert event3 == event2

    # just running this used to trigger a bug
    list(m(kin, 1))


@pytest.mark.filterwarnings("error")
def test_warning_when_creating_two_instances():
    m1 = Pythia6()

    with pytest.raises(RuntimeWarning):
        Pythia6()

    del m1

    # no warning
    Pythia6()
