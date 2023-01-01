from impy.kinematics import CenterOfMass
from impy.models import EposLHC
from impy.constants import GeV
import numpy as np
from numpy.testing import assert_allclose
from .util import reference_charge
import pytest
from particle import literals as lp
from functools import lru_cache


def run_pp_collision():
    m = EposLHC(seed=4)
    m.set_stable("pi0", True)
    kin = CenterOfMass(10 * GeV, "proton", "proton")
    for event in m(kin, 1):
        pass
    return event


def run_ab_collision():
    kin = CenterOfMass(10 * GeV, (4, 2), (12, 6))
    m = EposLHC(seed=1)
    for event in m(kin, 1):
        pass
    return event


@pytest.fixture
@lru_cache(maxsize=1)
def event():
    # must cache a copy and not a view
    return run_pp_collision().copy()


@pytest.fixture
@lru_cache(maxsize=1)
def event_ion():
    # must cache a copy and not a view
    return run_ab_collision().copy()


def test_impact_parameter(event, event_ion):
    # EPOS defines an impact parameter also for pp
    assert event.impact_parameter > 0
    assert event_ion.impact_parameter > 0


def test_n_wounded(event):
    # TODO Pythia8 returns (0, 0) for pp collision, perhaps unify the response
    assert event.n_wounded == (1, 1)


@pytest.mark.xfail(reason="FIXME: n_wounded always seems to return (1, 1) for EPOS")
def test_n_wounded_ion(event_ion):
    assert event_ion.n_wounded[0] > 1
    assert event_ion.n_wounded[1] > 1


def run_cross_section(p1, p2):
    m = EposLHC(seed=1)
    kin = CenterOfMass(10 * GeV, p1, p2)
    return m.cross_section(kin)


def test_cross_section():
    c = run_cross_section("p", "p")
    assert_allclose(c.total, 38.2, atol=0.1)
    assert_allclose(c.inelastic, 30.7, atol=0.1)
    assert_allclose(c.elastic, 7.4, atol=0.1)
    assert_allclose(c.diffractive_xb, 1.6, atol=0.1)
    assert_allclose(c.diffractive_ax, 1.6, atol=0.1)
    assert_allclose(c.diffractive_xx, 9.9, atol=0.1)
    assert c.diffractive_axb == 0
    assert_allclose(
        c.non_diffractive,
        c.inelastic - c.diffractive_xb - c.diffractive_ax - c.diffractive_xx,
    )


def test_charge(event):
    expected = reference_charge(event.pid)
    # skip internal particles unknown to reference_charge
    ma = np.isnan(expected)
    # EPOS has lots of unknown generator-specific particles.
    # For these, NaN is returned. We just check here whether
    # the fraction is not 100 %, that cannot be, since at least
    # final state particles cannot be generator-specific.
    assert np.mean(ma) < 0.8
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


@pytest.mark.parametrize("stable", (False, True))
def test_set_stable(stable):
    m = EposLHC(seed=4)
    pid = lp.pi_0.pdgid

    m.set_stable("pi_0", stable)

    kin = CenterOfMass(10 * GeV, "proton", "proton")
    for event in m(kin, 1):
        pass

    final = event.final_state()
    if stable:
        assert np.sum(final.pid == pid) > 0
    else:
        assert np.sum(final.pid == pid) == 0

    del m

    # check that defaults are restored for new instance
    m = EposLHC(seed=4)
    kin = CenterOfMass(10 * GeV, "proton", "proton")
    for event in m(kin, 1):
        pass

    final = event.final_state()
    assert np.sum(final.pid == pid) == 0
