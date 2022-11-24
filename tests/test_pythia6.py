# These tests check Pythia6 but also the MCEvent in general.
# It is not necessary to duplicate all tests for every model.
from impy.kinematics import CenterOfMass
from impy.models import Pythia6
from impy.constants import GeV, TeV
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from .util import reference_charge, run_in_separate_process
import pytest
import pickle
from particle import literals as lp
from functools import lru_cache


def test_name():
    assert Pythia6.name == "Pythia"
    assert Pythia6.version == "6.428"
    assert Pythia6.label == "Pythia-6.428"


def run_name():
    evt_kin = CenterOfMass(1 * TeV, "p", "p")
    m = Pythia6(evt_kin, seed=4)
    assert m.label == "Pythia-6.428"


def test_instance_name():
    run_in_separate_process(run_name)


def run_pp_collision():
    evt_kin = CenterOfMass(1 * TeV, "p", "p")
    m = Pythia6(evt_kin, seed=4)
    m.set_stable(lp.pi_0.pdgid, False)  # needed to get nonzero vertices
    for event in m(1):
        pass
    return event  # MCEvent is restored as EventData


def run_cross_section(p1, p2):
    evt_kin = CenterOfMass(10 * GeV, p1, p2)
    m = Pythia6(evt_kin, seed=1)
    return m.cross_section()


@pytest.fixture
@lru_cache(maxsize=1)
def event():
    return run_in_separate_process(run_pp_collision)


def test_cross_section():
    c = run_in_separate_process(run_cross_section, "p", "p")
    assert_allclose(c.total, 38.4, atol=0.1)
    assert_allclose(c.inelastic, 31.4, atol=0.1)
    assert_allclose(c.elastic, 7.0, atol=0.1)
    assert_allclose(c.diffractive_xb, 2.6, atol=0.1)
    assert_allclose(c.diffractive_ax, 2.6, atol=0.1)
    assert_allclose(c.diffractive_xx, 0.9, atol=1)
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


def run_is_view():
    evt_kin = CenterOfMass(10 * GeV, 2212, 2212)

    m = Pythia6(evt_kin, seed=1)
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


def run_pickle():
    evt_kin = CenterOfMass(10 * GeV, 2212, 2212)

    m = Pythia6(evt_kin, seed=1)
    for event in m(1):
        pass

    s = pickle.dumps(event)
    event2 = pickle.loads(s)

    assert event == event2


def test_pickle(event):
    run_in_separate_process(run_pickle)


def run_pp_collision_copy():
    from impy.models.pythia6 import PYTHIA6Event
    from impy.common import EventData

    evt_kin = CenterOfMass(1 * TeV, 2212, 2212)
    m = Pythia6(evt_kin, seed=4)
    m.set_stable(lp.pi_0.pdgid, False)  # needed to get nonzero vertices
    for event in m(1):
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
    list(m(1))


def test_event_copy():
    run_in_separate_process(run_pp_collision_copy)
