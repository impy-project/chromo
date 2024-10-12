import numpy as np
from numpy.testing import assert_allclose

import chromo
from chromo.constants import GeV
from chromo.util import naneq

from .util import run_in_separate_process


def run_cross_section(p1, p2):
    evt_kin = chromo.kinematics.CenterOfMass(10 * GeV, p1, p2)
    m = chromo.models.DpmjetIII193(evt_kin, seed=1)
    return m.cross_section(max_info=True)


def run_cross_section_ntrials():
    event_kin = chromo.kinematics.FixedTarget(1e4, "proton", "O16")
    event_generator = chromo.models.DpmjetIII193(event_kin)

    air = chromo.util.CompositeTarget([("N", 0.78), ("O", 0.22)])
    event_kin = chromo.kinematics.FixedTarget(1e3, "proton", air)

    default_precision = 1000
    # Check the default precision
    assert event_generator.glauber_trials == default_precision

    # Set a new one
    other_precision = 58
    event_generator.glauber_trials = other_precision
    assert event_generator.glauber_trials == other_precision

    trials = 10

    # With small precision
    event_generator.glauber_trials = 1
    cross_section_run1 = np.empty(trials, dtype=np.float64)
    for i in range(trials):
        cross_section_run1[i] = event_generator.cross_section(
            event_kin, max_info=True
        ).prod

    # With default precision
    event_generator.glauber_trials = 1000
    cross_section_run2 = np.empty(trials, dtype=np.float64)
    for i in range(trials):
        cross_section_run2[i] = event_generator.cross_section(
            event_kin, max_info=True
        ).prod

    # Standard deviation should be large for small precision
    assert np.std(cross_section_run1) > np.std(cross_section_run2)


def test_cross_section_ntrials():
    run_in_separate_process(run_cross_section_ntrials)


def test_cross_section_pp():
    c = run_in_separate_process(run_cross_section, "p", "p")
    # These are the expected rounded numbers from the DPMJET
    # for pp at 10 GeV
    assert_allclose(c.total, 38.9, atol=0.1)
    assert_allclose(c.inelastic, 32.1, atol=0.1)
    assert_allclose(c.elastic, 6.9, atol=0.1)
    assert_allclose(
        c.non_diffractive,
        c.inelastic,
    )
    naneq(c.diffractive_ax, np.nan)
    naneq(c.diffractive_xb, np.nan)
    naneq(c.diffractive_xx, np.nan)
    naneq(c.diffractive_axb, np.nan)


def test_cross_section_pA():
    c = run_in_separate_process(run_cross_section, "p", "O16")
    assert_allclose(c.total, 446.1, atol=0.1)
    assert_allclose(c.inelastic, 328.0, atol=0.1)
    assert_allclose(c.elastic, 118.1, atol=0.1)
    assert_allclose(c.prod, 298.7, atol=0.1)
    assert_allclose(c.quasielastic, 144.6, atol=0.1)
    assert_allclose(
        c.non_diffractive,
        c.inelastic,
    )
    naneq(c.diffractive_ax, np.nan)
    naneq(c.diffractive_xb, np.nan)
    naneq(c.diffractive_xx, np.nan)
    naneq(c.diffractive_axb, np.nan)
