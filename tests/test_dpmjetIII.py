import sys

import numpy as np
import pytest
from numpy.testing import assert_allclose

import chromo
from chromo.constants import GeV
from chromo.util import get_all_models, naneq

from .util import run_in_separate_process

pytestmark = pytest.mark.skipif(
    sys.platform == "win32", reason="DPMJETIII19x not build on windows"
)


def get_dpmjets(no307=False):
    """Get the list of all DPMJets"""
    return [
        cl
        for cl in get_all_models()
        if cl.pyname.startswith("Dpmjet") and (not no307 or cl.pyname != "DpmjetIII307")
    ]


def run_cross_section(p1, p2, model):
    evt_kin = chromo.kinematics.CenterOfMass(10 * GeV, p1, p2)
    m = model(evt_kin, seed=1)
    return m.cross_section(max_info=True)


def run_three_events(p1, model):
    chromo.debug_level = 1
    evt_kin = chromo.kinematics.CenterOfMass(100 * GeV, p1, "O16")
    m = model(evt_kin, seed=1)
    for evt in m(3):
        evt = evt.final_state()
        assert len(evt.en) > 0


def run_init_energy_guard(model):
    # DPMJET tabulates hadron-nucleon cross sections in dt_init only up
    # to the initialization energy; requesting higher-energy kinematics
    # afterwards must raise instead of silently extrapolating.
    evt_kin = chromo.kinematics.FixedTarget(1e3, "proton", "O16")
    m = model(evt_kin, seed=1)
    # at or below the initialization energy: allowed
    m.cross_section(chromo.kinematics.FixedTarget(1e2, "proton", "O16"))
    try:
        m.cross_section(chromo.kinematics.FixedTarget(1e4, "proton", "O16"))
    except ValueError:
        return True
    return False


@pytest.mark.parametrize("model", get_dpmjets())
def test_init_energy_guard(model):
    assert run_in_separate_process(run_init_energy_guard, model)


def run_cross_section_ntrials(model):
    event_kin = chromo.kinematics.FixedTarget(1e4, "proton", "O16")
    event_generator = model(event_kin)

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


@pytest.mark.parametrize("model", get_dpmjets())
def test_cross_section_ntrials(model):
    run_in_separate_process(run_cross_section_ntrials, model)


@pytest.mark.parametrize("model", get_dpmjets(no307=True))
def test_cross_section_pp(model):
    c = run_in_separate_process(run_cross_section, "p", "p", model)
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


@pytest.mark.parametrize("model", get_dpmjets(no307=True))
def test_cross_section_pA(model):
    c = run_in_separate_process(run_cross_section, "p", "O16", model)
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


def run_cross_section_pA_energy_dependence(model):
    # Issue #242: the production cross section for h+A must be calculated
    # for the queried kinematics and not be the stale value left in the
    # tables by the initialization at the highest energy.
    gen = model(chromo.kinematics.FixedTarget(1e5, "proton", "O16"), seed=1)

    prod = {
        plab: gen.cross_section(
            chromo.kinematics.FixedTarget(plab, "proton", "O16")
        ).prod
        for plab in (1e2, 1e5)
    }

    # physical and energy dependent
    assert all(np.isfinite(p) and p > 0 for p in prod.values()), str(prod)
    assert prod[1e2] < prod[1e5]

    # consistent with the full Glauber MC estimate at the same energies
    for plab in (1e2, 1e5):
        full = gen.cross_section(
            chromo.kinematics.FixedTarget(plab, "proton", "O16"), max_info=True
        )
        assert_allclose(prod[plab], full.prod, rtol=0.05)


@pytest.mark.parametrize("model", get_dpmjets())
def test_cross_section_pA_energy_dependence(model):
    run_in_separate_process(run_cross_section_pA_energy_dependence, model)


def run_prod_cs_event_stream(model, prod_queries):
    gen = model(chromo.kinematics.FixedTarget(1e5, "proton", "O16"), seed=1)
    if prod_queries:
        for plab in (1e2, 1e4):
            gen.cross_section(chromo.kinematics.FixedTarget(plab, "proton", "O16"))
    return [(len(evt.final_state()), np.sum(evt.final_state().en)) for evt in gen(3)]


@pytest.mark.parametrize("model", get_dpmjets())
def test_cross_section_pA_rng_neutral(model):
    # prod-only Glauber runs for cross-section queries must save/restore
    # the RNG state, so the event generation stream is unaffected
    ref = run_in_separate_process(run_prod_cs_event_stream, model, False)
    with_queries = run_in_separate_process(run_prod_cs_event_stream, model, True)
    for (n1, e1), (n2, e2) in zip(ref, with_queries):
        assert n1 == n2
        assert_allclose(e1, e2, rtol=1e-10)


def get_model_projectile_combinations():
    """Get combinations of DPMJET models and their non-nuclei projectiles with PDG ID < 6000"""
    return [
        (model, int(pid))
        for model in get_dpmjets()
        for pid in getattr(model.projectiles, "_other", set())
        if int(pid) < 6000
    ]


@pytest.mark.parametrize("model,p1", get_model_projectile_combinations())
def test_projectile_list(model, p1):
    run_in_separate_process(run_three_events, p1, model)
