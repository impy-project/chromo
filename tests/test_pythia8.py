from chromo.kinematics import CenterOfMass
from chromo.models import Pythia8
from chromo.constants import GeV, long_lived
from chromo.util import name2pdg
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from .util import reference_charge
import pytest
from functools import lru_cache
import sys


pytestmark = pytest.mark.skipif(
    sys.platform == "win32", reason="Pythia8 does not run on windows"
)


def run_collision(energy, p1, p2):
    evt_kin = CenterOfMass(energy, p1, p2)
    m = Pythia8(evt_kin, seed=4)
    for event in m(1):
        pass
    return event


def run_cross_section(energy, p1, p2):
    evt_kin = CenterOfMass(energy, p1, p2)
    m = Pythia8(evt_kin, seed=1)
    return m.cross_section()


@pytest.fixture
@lru_cache(maxsize=1)  # Pythia8 initialization is very slow
def event():
    return run_collision(10 * GeV, "p", "p")


def test_impact_parameter(event):
    assert np.isnan(event.impact_parameter)


def test_n_wounded(event):
    # TODO EPOS returns (1, 1) for pp collision, perhaps unify the response
    assert event.n_wounded == (0, 0)


def test_cross_section():
    c = run_cross_section(10 * GeV, "p", "p")
    assert_allclose(c.total, 38.4, atol=0.1)
    assert_allclose(c.inelastic, 31.3, atol=0.1)
    assert_allclose(c.elastic, 7.1, atol=0.1)
    assert_allclose(c.diffractive_xb, 2.6, atol=0.1)
    assert_allclose(c.diffractive_ax, 2.6, atol=0.1)
    assert_allclose(c.diffractive_xx, 0.8, atol=0.1)
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


def test_daughters(event):
    assert event.daughters.shape == (len(event), 2)
    # some particles have no daughters
    assert sum(x[0] == -1 and x[1] == -1 for x in event.daughters) > 0

    # somes particles have single daughters (elastic scattering in parton shower)
    assert sum(x[0] > -1 and x[1] == -1 for x in event.daughters) > 0

    # some particles have multiple daughters
    assert sum(x[0] > -1 and x[1] > -1 for x in event.daughters) > 0


def test_mothers(event):
    n = len(event)
    assert event.mothers.shape == (n, 2)
    # same particles have no mothers
    assert sum(x[0] == -1 and x[1] == -1 for x in event.mothers) > 0

    # most particles have a single mother
    assert sum(x[0] > -1 and x[1] == -1 for x in event.mothers) > 0

    # some particles have multiple mothers
    assert sum(x[0] > -1 and x[1] > -1 for x in event.mothers) > 0

    assert sum(x[1] >= -1 and x[1] < len(event) for x in event.mothers) == n


@pytest.mark.parametrize("projectile", ("pi-", "K+", "p", "He"))
def test_nuclear_collision(projectile):
    # This test is very slow the first time when Pythia is computing
    # some cached results. Afterwards it runs super fast.
    ecm = 2000 * GeV
    event = run_collision(ecm, projectile, "He")
    assert event.pid[0] == name2pdg(projectile)
    assert event.pid[1] == 1000020040
    if projectile == "p":
        assert_allclose(event.en[0], ecm / 2)
        assert_allclose(event.en[1], ecm * 2, rtol=1e-2)
    elif projectile == "He":
        assert_allclose(event.en[0], ecm * 2, rtol=1e-2)
        assert_allclose(event.en[1], ecm * 2, rtol=1e-2)
    assert event.impact_parameter > 0
    if projectile == "He":
        assert event.n_wounded[0] >= 1
    else:
        assert event.n_wounded[0] == 1
    assert event.n_wounded[1] > 0
    apid = np.abs(event.final_state_charged().pid)
    assert np.sum(apid == 211) > 10


def test_photo_hadron_collision():
    event = run_collision(100 * GeV, "gamma", "p")
    event.pid[0] == 22
    event.pid[1] == 2212
    apid = np.abs(event.final_state_charged().pid)
    assert np.sum(apid == 211) > 0


def test_changing_beams_proton():
    evt_kin = CenterOfMass(10 * GeV, "p", "p")
    m = Pythia8(evt_kin, seed=1)
    for event in m(1):
        assert_allclose(event.en[:2], 5 * GeV)
    m.kinematics = CenterOfMass(100 * GeV, "p", "p")
    for event in m(1):
        assert_allclose(event.en[:2], 50 * GeV)


def test_event(event):
    assert_equal(event.pid[:2], (2212, 2212))


def test_elastic():
    evt_kin = CenterOfMass(10 * GeV, "p", "p")
    m = Pythia8(evt_kin, seed=1, config=["SoftQCD:elastic=on"])
    for event in m(10):
        assert len(event) == 4
        assert_equal(event.pid, [2212] * 4)


def test_gamma_p():
    evt_kin = CenterOfMass(10 * GeV, "gamma", "p")
    m = Pythia8(evt_kin, seed=1)
    for event in m(10):
        assert len(event) > 2


@pytest.mark.parametrize("seed", (None, 0, 1, int(1e10)))
def test_seed(seed):
    evt_kin = CenterOfMass(10 * GeV, "p", "p")
    m = Pythia8(evt_kin, seed=seed)
    if seed is None:
        assert m.seed >= 0
    else:
        assert m.seed == seed


def test_get_stable():
    evt_kin = CenterOfMass(10 * GeV, "p", "p")
    m = Pythia8(evt_kin, seed=1)
    assert m._get_stable() == set(long_lived)


def test_gp():
    evt = run_collision(100 * GeV, "gamma", "p")
    assert len(evt) > 2


@pytest.mark.xfail(reason="no rest frame for gamma-gamma")
def test_gg():
    evt = run_collision(100 * GeV, "gamma", "gamma")
    assert len(evt) > 2
