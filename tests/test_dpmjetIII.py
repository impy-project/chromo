import sys

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

import chromo
from chromo.constants import GeV
from chromo.util import get_all_models, naneq

from .util import baryon_number as _baryon_number
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
    return [(len(evt.final_state()), np.sum(evt.final_state().en)) for evt in gen(6)]


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


def run_photon_on_nucleus(model, target, elab):
    chromo.debug_level = 1
    evt_kin = chromo.kinematics.FixedTarget(elab * GeV, "gamma", target)
    m = model(evt_kin, seed=1)
    assert 22 in model.projectiles
    xs = m.cross_section()
    assert 0 < xs.prod < 100, f"photon production xs out of range: {xs.prod}"
    for evt in m(3):
        assert evt.pid[0] == 22
        assert len(evt.final_state().en) > 0
    # max_info runs the Glauber MC and consumes the Fortran RNG state,
    # so it must come last and forbid further event generation
    xs = m.cross_section(max_info=True)
    assert 0 < xs.total < 100
    assert xs.elastic < xs.total
    assert np.isclose(xs.inelastic, xs.total - xs.elastic)
    try:
        next(iter(m(1)))
    except RuntimeError:
        pass
    else:
        return False
    return True


@pytest.mark.parametrize("model", get_dpmjets())
@pytest.mark.parametrize("target", ["O16", "Fe56"])
def test_photon_nucleus(model, target):
    assert run_in_separate_process(run_photon_on_nucleus, model, target, 1e4)


def run_photon_on_nucleon(model, target):
    m = model(chromo.kinematics.FixedTarget(1e4 * GeV, "gamma", target), seed=1)
    xs = m.cross_section()
    for evt in m(3):
        assert evt.pid[0] == 22
        assert len(evt.final_state().en) > 0
    return xs.inelastic


@pytest.mark.parametrize("target", ["p", "n"])
def test_dpmjet193_photon_nucleon(target):
    from chromo.models import DpmjetIII193

    sine = run_in_separate_process(run_photon_on_nucleon, DpmjetIII193, target)
    # PHOJET sigma_inel(gamma p) at sqrt(s) = 137 GeV
    assert_allclose(sine, 0.1516, rtol=0.01)


def run_photon_on_nucleon_rejected(model, target):
    try:
        model(chromo.kinematics.FixedTarget(1e3 * GeV, "gamma", target), seed=1)
    except ValueError:
        return True
    return False


@pytest.mark.parametrize("target", ["p", "n"])
def test_dpmjet307_photon_nucleon_rejected(target):
    from chromo.models import DpmjetIII307

    assert run_in_separate_process(run_photon_on_nucleon_rejected, DpmjetIII307, target)


def run_first_event(model, kin, randomize):
    gen = model(kin, seed=1)
    gen.randomize_azimuth = randomize
    return next(gen(1)).copy()


@pytest.mark.parametrize("model", get_dpmjets())
def test_randomize_azimuth_is_rigid_rotation(model):
    # same seed: first events differ only by a rotation around the beam axis
    kin = chromo.kinematics.CenterOfMass(100 * GeV, "O", "O16")
    rot = run_in_separate_process(run_first_event, model, kin, True)
    ref = run_in_separate_process(run_first_event, model, kin, False)
    assert_equal(rot.pid, ref.pid)
    assert_equal(rot.status, ref.status)
    for attr in ("en", "pz", "m", "vz"):
        assert_allclose(getattr(rot, attr), getattr(ref, attr), rtol=1e-8)
    pt_ref = np.hypot(ref.px, ref.py)
    assert_allclose(np.hypot(rot.px, rot.py), pt_ref, rtol=1e-8)
    assert_allclose(np.hypot(rot.vx, rot.vy), np.hypot(ref.vx, ref.vy), rtol=1e-8)
    dphi = np.arctan2(rot.py, rot.px) - np.arctan2(ref.py, ref.px)
    dphi = (dphi + np.pi) % (2 * np.pi) - np.pi
    sel = pt_ref > 1e-6
    assert_allclose(dphi[sel], dphi[sel][0], atol=1e-6)
    assert np.max(np.abs(rot.px[sel] - ref.px[sel])) > 1e-6


@pytest.mark.parametrize("model", get_dpmjets())
def test_randomize_azimuth_on_by_default(model):
    assert model.randomize_azimuth is True


@pytest.mark.parametrize("model", get_dpmjets())
def test_randomize_azimuth_skips_hadron_nucleon(model):
    kin = chromo.kinematics.CenterOfMass(100 * GeV, "p", "p")
    a = run_in_separate_process(run_first_event, model, kin, True)
    b = run_in_separate_process(run_first_event, model, kin, False)
    assert_equal(a.pid, b.pid)
    assert_allclose(a.en, b.en, rtol=1e-8)
    assert_allclose(a.px, b.px, rtol=1e-8)


def run_remnant_conservation(model, p1, p2):
    import numpy as np

    from chromo.kinematics import FixedTarget

    from .util import charge_number

    kin = FixedTarget(1e4, p1, p2)
    a = (kin.p1.A or 0) + (kin.p2.A or 0)
    z = (kin.p1.Z or 0) + (kin.p2.Z or 0)
    generator = model(kin, seed=1)
    for event in generator(3):
        assert not np.any((event.status == 5) & (event.daughters[:, 0] != -1))
        frags = event.final_state_with_nucl_frag()
        assert sum(_baryon_number(pid) for pid in frags.pid) == a
        assert sum(charge_number(pid) for pid in frags.pid) == z


@pytest.mark.parametrize("model", get_dpmjets(no307=False))
@pytest.mark.parametrize("p1,p2", [("p", "O16"), ("p", "Pb208"), ("O16", "O16")])
def test_remnant_conservation(model, p1, p2):
    run_in_separate_process(run_remnant_conservation, model, p1, p2)
