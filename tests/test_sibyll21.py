from chromo.kinematics import CenterOfMass
from chromo.models import Sibyll21
from chromo.constants import TeV, GeV
from numpy.testing import assert_allclose, assert_equal
import numpy as np
from .util import reference_charge
import pytest
from particle import Particle


def event_run():
    m = Sibyll21(seed=1)
    kin = CenterOfMass(10 * TeV, "p", "p")
    for event in m(kin, 1):
        pass
    return event.copy()  # copy is pickleable


@pytest.fixture
def event():
    return event_run()


@pytest.mark.xfail(reason="needs a fix to sibylls ICHP table")
def test_charge(event):
    expected = reference_charge(event.pid)
    assert_allclose(event.charge, expected)


def test_children(event):
    assert event.children is None


def test_parents(event):
    # check that there are particles with a single parent
    # and that parent is short-lived or an initial particle

    d = event.parents[:, 1] - event.parents[:, 0]
    ma = (d == 0) & (event.parents[:, 0] > 0)

    assert np.sum(ma) > 0

    idx = event.parents[ma, 0] - 1

    # remove beam particles (if any)
    ma = event.status[idx] != 3
    idx = idx[ma]

    # remaining particles should be short-lived or
    # bookkeeping particles like K0
    for pid in event.pid[idx]:
        p = Particle.from_pdgid(pid)
        assert p.ctau is None or p.ctau < 1


def test_vertex(event):
    # no vertex info in Sibyll21
    assert_equal(event.vx, 0)
    assert_equal(event.vy, 0)
    assert_equal(event.vz, 0)
    assert_equal(event.vt, 0)


def test_cross_section():
    m = Sibyll21()
    kin = CenterOfMass(10 * GeV, "p", "p")
    c = m.cross_section(kin)
    assert_allclose(c.total, 38.4, atol=0.1)
    assert_allclose(c.inelastic, 30.9, atol=0.1)
    assert_allclose(c.elastic, 7.4, atol=0.1)
    assert_allclose(c.diffractive_xb, 2.9, atol=0.1)
    assert_allclose(c.diffractive_ax, 2.9, atol=0.1)
    assert_allclose(c.diffractive_xx, 0.8, atol=0.1)
    assert c.diffractive_axb == 0
    assert_allclose(
        c.non_diffractive,
        c.inelastic - c.diffractive_xb - c.diffractive_ax - c.diffractive_xx,
    )
