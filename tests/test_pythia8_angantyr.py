import sys
from functools import lru_cache

import numpy as np
import pytest
from numpy.testing import assert_allclose

from chromo.constants import GeV
from chromo.kinematics import CenterOfMass, FixedTarget
from chromo.models import Pythia8Angantyr

from .util import reference_charge, run_in_separate_process

pytestmark = pytest.mark.skipif(
    sys.platform == "win32", reason="Pythia8 does not run on windows"
)


def run_angantyr_collision(energy, p1, p2, fixed_target=False):
    if fixed_target:
        kin = FixedTarget(energy, p1, p2)
    else:
        kin = CenterOfMass(energy, p1, p2)
    gen = Pythia8Angantyr(kin, seed=1)
    for event in gen(1):
        pass
    return event


def run_angantyr_cross_section(energy, p1, p2):
    kin = CenterOfMass(energy, p1, p2)
    gen = Pythia8Angantyr(kin, seed=1)
    # Glauber cross sections accumulate during event generation
    for _ in gen(100):
        pass
    return gen.cross_section()


def run_angantyr_changing_kinematics():
    """Run two different energies with the same generator instance."""
    kin1 = CenterOfMass(100, "p", "O16")
    gen = Pythia8Angantyr(kin1, seed=1)
    for event in gen(1):
        pass
    # Change kinematics — should use setBeamIDs + setKinematics, not re-init
    gen.kinematics = CenterOfMass(200, "p", "O16")
    for event in gen(1):
        pass
    return event


@pytest.fixture
@lru_cache(maxsize=1)
def event_pp():
    return run_in_separate_process(run_angantyr_collision, 100, "p", "p")


@pytest.fixture
@lru_cache(maxsize=1)
def event_p_O16():
    return run_in_separate_process(run_angantyr_collision, 100, "p", "O16")


def test_pp_collision(event_pp):
    fs = event_pp.final_state_charged()
    assert len(fs) > 0, "No final-state charged particles in p+p"


def test_p_O16_collision(event_p_O16):
    fs = event_p_O16.final_state_charged()
    apid = np.abs(fs.pid)
    assert np.sum(apid == 211) > 0, "No charged pions in final state"
    assert len(fs) > 2, "Too few final-state particles"


def test_p_N14_cross_section():
    cs = run_in_separate_process(run_angantyr_cross_section, 100, "p", "N14")
    # p+N14 inelastic at 100 GeV CMS should be ~300 mb
    assert cs.inelastic > 100, f"Cross section too small: {cs.inelastic} mb"
    assert cs.inelastic < 2000, f"Cross section too large: {cs.inelastic} mb"


def test_impact_parameter(event_p_O16):
    assert np.isfinite(event_p_O16.impact_parameter), (
        f"Impact parameter not finite: {event_p_O16.impact_parameter}"
    )


def test_n_wounded(event_p_O16):
    nw = event_p_O16.n_wounded
    assert nw[0] >= 1 or nw[1] >= 1, f"No wounded nucleons: {nw}"


def test_charge(event_p_O16):
    expected = reference_charge(event_p_O16.pid)
    ma = np.isnan(expected)
    assert np.mean(ma) < 0.1
    event_p_O16.charge[ma] = np.nan
    assert_allclose(event_p_O16.charge, expected)


def test_daughters(event_p_O16):
    assert event_p_O16.daughters.shape == (len(event_p_O16), 2)


def test_mothers(event_p_O16):
    assert event_p_O16.mothers.shape == (len(event_p_O16), 2)


def test_changing_kinematics():
    event = run_in_separate_process(run_angantyr_changing_kinematics)
    fs = event.final_state_charged()
    assert len(fs) > 0


def test_fixed_target_mode():
    event = run_in_separate_process(
        run_angantyr_collision, 100 * GeV, "p", "O16", True
    )
    fs = event.final_state_charged()
    assert len(fs) > 2
