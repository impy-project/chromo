import sys
from functools import lru_cache

import numpy as np
import pytest
from numpy.testing import assert_allclose

from chromo.constants import GeV
from chromo.kinematics import FixedTarget
from chromo.models import Pythia8Cascade

from .util import reference_charge, run_in_separate_process

pytestmark = pytest.mark.skipif(
    sys.platform == "win32", reason="Pythia8 does not run on windows"
)


def run_cascade_collision(energy, p1, p2):
    kin = FixedTarget(energy, p1, p2)
    gen = Pythia8Cascade(kin, seed=1)
    for event in gen(1):
        pass
    return event


def run_cascade_cross_section(energy, p1, p2):
    kin = FixedTarget(energy, p1, p2)
    gen = Pythia8Cascade(kin, seed=1)
    return gen.cross_section()


@pytest.fixture
@lru_cache(maxsize=1)
def event_p_O16():
    return run_in_separate_process(run_cascade_collision, 100 * GeV, "p", "O16")


@pytest.fixture
@lru_cache(maxsize=1)
def event_He_O16():
    return run_in_separate_process(run_cascade_collision, 100 * GeV, "He", "O16")


def test_cross_section_p_N14():
    cs = run_in_separate_process(run_cascade_cross_section, 100 * GeV, "p", "N14")
    # p+N14 inelastic cross section at ~100 GeV/c should be several hundred mb
    assert cs.inelastic > 100, f"Cross section too small: {cs.inelastic} mb"
    assert cs.inelastic < 2000, f"Cross section too large: {cs.inelastic} mb"


def test_cross_section_p_O16():
    cs = run_in_separate_process(run_cascade_cross_section, 100 * GeV, "p", "O16")
    # p+O16 inelastic cross section at ~100 GeV/c should be several hundred mb
    assert cs.inelastic > 100, f"Cross section too small: {cs.inelastic} mb"
    assert cs.inelastic < 2000, f"Cross section too large: {cs.inelastic} mb"


def test_cross_section_p_Pb208():
    cs = run_in_separate_process(run_cascade_cross_section, 100 * GeV, "p", "Pb208")
    # p+Pb inelastic cross section should be ~3000 mb
    assert cs.inelastic > 500, f"Cross section too small: {cs.inelastic} mb"
    assert cs.inelastic < 8000, f"Cross section too large: {cs.inelastic} mb"


def test_n_wounded_p_O16(event_p_O16):
    # At least one target nucleon should be wounded
    assert event_p_O16.n_wounded[0] >= 1


def test_final_state_p_O16(event_p_O16):
    fs = event_p_O16.final_state_charged()
    apid = np.abs(fs.pid)
    assert np.sum(apid == 211) > 0, "No charged pions in final state"
    assert len(fs) > 2, "Too few final-state particles"


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


def test_nuclear_projectile_He_O16(event_He_O16):
    fs = event_He_O16.final_state_charged()
    assert len(fs) > 2


def test_nuclear_projectile_Fe_O16():
    event = run_in_separate_process(run_cascade_collision, 100 * GeV, "Fe56", "O16")
    fs = event.final_state_charged()
    assert len(fs) > 2
