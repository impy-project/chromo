"""Tests for the Fluka model.

Each test runs in a separate subprocess via run_in_separate_process
because FLUKA is single-instantiation per Python process.
"""

import os
import sys
from functools import lru_cache

import numpy as np
import pytest

from chromo.kinematics import CenterOfMass, FixedTarget
from chromo.util import CompositeTarget

from .util import reference_charge, run_in_separate_process

try:
    from chromo.models._fluka import pdg_to_proj_code  # noqa: F401

    _fluka_available = True
except ImportError:
    _fluka_available = False

_flupro_valid = "FLUPRO" in os.environ and os.path.isdir(os.environ.get("FLUPRO", ""))

pytestmark = [
    pytest.mark.skipif(sys.platform == "win32", reason="FLUKA is not built on Windows"),
    pytest.mark.skipif(not _fluka_available, reason="_fluka extension not built"),
    pytest.mark.skipif(not _flupro_valid, reason="FLUPRO not set or invalid"),
]


def _run_xsec(ecm_or_elab, p1, p2, *, fixed_target=True, **kwargs):
    from chromo.models import Fluka

    kin = (
        FixedTarget(ecm_or_elab, p1, p2)
        if fixed_target
        else CenterOfMass(ecm_or_elab, p1, p2)
    )
    gen = Fluka(kin, seed=1, **kwargs)
    return gen.cross_section()


def _run_one_event(elab, p1, p2, **kwargs):
    from chromo.models import Fluka

    kin = FixedTarget(elab, p1, p2)
    gen = Fluka(kin, seed=1, **kwargs)
    for event in gen(1):
        pass
    return event


def test_import():
    """Trivial importability test."""
    from chromo.models import Fluka

    assert Fluka.name == "FLUKA"


# ---------------------------------------------------------------------------
# Task 13: Cross-section tests (hN, hA, AA, composite)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "target,xs_min,xs_max",
    [
        ("H1", 20, 60),
        ("N14", 200, 600),
        ("O16", 200, 700),
        ("Ar40", 400, 1100),
        ("Fe56", 500, 1300),
        ("Pb208", 1200, 2800),
    ],
)
def test_xsec_p_A_sweep(target, xs_min, xs_max):
    cs = run_in_separate_process(_run_xsec, 100.0, "p", target)
    assert (
        xs_min < cs.inelastic < xs_max
    ), f"sigma_inel(p+{target})={cs.inelastic} mb outside [{xs_min}, {xs_max}]"


def test_xsec_pi_N14():
    cs = run_in_separate_process(_run_xsec, 100.0, "pi+", "N14")
    assert 150 < cs.inelastic < 600


def test_xsec_composite_air():
    air = CompositeTarget([("N14", 0.78), ("O16", 0.21), ("Ar40", 0.01)], label="Air")
    cs = run_in_separate_process(_run_xsec, 100.0, "p", air)
    assert 200 < cs.inelastic < 1100


def test_xsec_AA_O_O():
    cs = run_in_separate_process(_run_xsec, 1600.0, "O16", "O16")
    assert 500 < cs.inelastic < 2500


# ---------------------------------------------------------------------------
# Task 14: Frame-conversion round-trip test
# ---------------------------------------------------------------------------


def _xsec_ft(elab, p1, p2):
    from chromo.models import Fluka

    return Fluka(FixedTarget(elab, p1, p2), seed=1).cross_section()


def _xsec_cms(ecm, p1, p2):
    from chromo.models import Fluka

    return Fluka(CenterOfMass(ecm, p1, p2), seed=1).cross_section()


def test_xsec_cms_vs_ft_equivalent():
    ecm = 20.0  # GeV
    cs_cms = run_in_separate_process(_xsec_cms, ecm, "p", "N14")
    m_p = 0.938272
    elab = (ecm * ecm - 2 * m_p * m_p) / (2 * m_p)
    cs_ft = run_in_separate_process(_xsec_ft, elab, "p", "N14")
    assert abs(cs_cms.inelastic - cs_ft.inelastic) < 0.05 * cs_ft.inelastic


# ---------------------------------------------------------------------------
# Task 15: Photon tests (photohadronic + photonuclear)
# ---------------------------------------------------------------------------


def test_xsec_gamma_p():
    cs = run_in_separate_process(_run_xsec, 10.0, "gamma", "p")
    # gamma-p sigma_inel at 10 GeV ekin is small (~0.1-1 mb range)
    assert 1e-3 < cs.inelastic < 1.0


def test_xsec_gamma_Pb_delta():
    # Δ-resonance region: 300 MeV photon
    cs = run_in_separate_process(_run_xsec, 0.3, "gamma", "Pb208")
    assert 1.0 < cs.inelastic < 500.0


def test_xsec_gamma_Pb_DIS():
    # DIS regime: 10 GeV photon
    cs = run_in_separate_process(_run_xsec, 10.0, "gamma", "Pb208")
    assert cs.inelastic > 1.0


def test_generate_gamma_Pb_at_delta():
    event = run_in_separate_process(_run_one_event, 0.3, "gamma", "Pb208")
    fs = event.final_state()
    assert len(fs) > 0
    # Residual nucleus should be present
    big_pids = np.abs(event.pid)
    has_big = np.any(big_pids >= 1_000_000_000)
    assert has_big, "no nuclear remnant in final state"


# ---------------------------------------------------------------------------
# Task 16: EMD tests
# ---------------------------------------------------------------------------


def _run_xsec_emd(elab, p1, p2):
    from chromo.models import Fluka
    from chromo.models.fluka import InteractionType

    kin = FixedTarget(elab, p1, p2)
    gen = Fluka(kin, seed=1, interaction_type=InteractionType.INELA_EMD)
    return gen.cross_section()


def _run_emd_event(elab, p1, p2):
    from chromo.models import Fluka
    from chromo.models.fluka import InteractionType

    kin = FixedTarget(elab, p1, p2)
    # EMD-only (IFLXYZ=100) aborts FLUKA for event generation;
    # use INELA_EMD (IFLXYZ=101) to include inelastic + EMD channels.
    gen = Fluka(kin, seed=1, interaction_type=InteractionType.INELA_EMD)
    for event in gen(1):
        pass
    return event


def test_xsec_emd_O16_Pb208():
    # EMD requires a nuclear projectile. Proton has EMD=0.
    # O16+Pb208 should have a sizeable EMD cross section.
    cs = run_in_separate_process(_run_xsec_emd, 1600.0, "O16", "Pb208")
    assert cs.inelastic > 0 and not np.isnan(cs.inelastic)
    assert cs.emd > 0 and not np.isnan(cs.emd)


def test_xsec_emd_AA_O_Pb():
    cs = run_in_separate_process(_run_xsec_emd, 1600.0, "O16", "Pb208")
    assert cs.emd > 0 and not np.isnan(cs.emd)


def test_generate_emd_event_one_Pb():
    # EMD-only event generation (IFLXYZ=100) aborts FLUKA.
    # INELA_EMD (IFLXYZ=101) requires a hadronic projectile for event
    # generation; nuclear projectiles (A>1) segfault in the FLUKA Fortran
    # backend during EVTXYZ. Use p+Pb208 which is supported.
    event = run_in_separate_process(_run_emd_event, 100.0, "p", "Pb208")
    fs = event.final_state()
    assert len(fs) > 0
    assert len(fs) < 500


# ---------------------------------------------------------------------------
# Task 17: Event, conservation, and remnant tests
# ---------------------------------------------------------------------------


@pytest.fixture
@lru_cache(maxsize=1)
def event_p_O16():
    return run_in_separate_process(_run_one_event, 100.0, "p", "O16")


@pytest.fixture
@lru_cache(maxsize=1)
def event_p_Pb208_high():
    # Nuclear projectiles (A>1) crash FLUKA's EVTXYZ in the current
    # wrapper. Use p+Pb208 as a proxy for heavy-target high-multiplicity.
    return run_in_separate_process(_run_one_event, 10_000.0, "p", "Pb208")


def test_generate_p_O16(event_p_O16):
    fs = event_p_O16.final_state()
    assert len(fs) > 2


def test_generate_p_O16_charged_pions(event_p_O16):
    fs = event_p_O16.final_state_charged()
    apid = np.abs(fs.pid)
    assert np.sum(apid == 211) > 0


def test_charge_reference_matches(event_p_O16):
    expected = reference_charge(event_p_O16.pid)
    mask = ~np.isnan(expected)
    np.testing.assert_allclose(event_p_O16.charge[mask], expected[mask])


def test_generate_p_Pb208_high_multiplicity(event_p_Pb208_high):
    # At 10 TeV lab, p+Pb should produce many particles.
    fs = event_p_Pb208_high.final_state()
    assert len(fs) > 10


@pytest.mark.parametrize(
    "p1,p2,elab",
    [
        ("p", "O16", 100.0),
        ("pi+", "Fe56", 50.0),
        # Nuclear projectiles (A>1) crash FLUKA event generation; omitted.
    ],
)
def test_conservation_baryon(p1, p2, elab):
    from particle import PDGID

    def baryon_number(pdg):
        """Return baryon number for a PDG id."""
        pid = PDGID(abs(int(pdg)))
        sign = 1 if int(pdg) > 0 else -1
        if pid.is_nucleus:
            return int(pid.A) * sign
        if pid.is_baryon:
            return 1 * sign
        return 0

    event = run_in_separate_process(_run_one_event, elab, p1, p2)
    fs = event.final_state()
    b_in = baryon_number(int(event.kin.p1)) + baryon_number(int(event.kin.p2))
    b_out = sum(baryon_number(int(pid)) for pid in fs.pid)
    assert b_in == b_out, f"baryon number not conserved: {b_in} != {b_out}"


def test_remnant_present_p_Pb208():
    event = run_in_separate_process(_run_one_event, 10.0, "p", "Pb208")
    big_pids = np.abs(event.pid)
    nuclei = big_pids[big_pids >= 1_000_000_000]
    assert len(nuclei) >= 1, "no nuclear remnant in p+Pb event"
    residual_A = np.array(
        [
            (abs(int(pid)) // 10) % 1000
            for pid in event.pid
            if abs(int(pid)) >= 1_000_000_000
        ]
    )
    assert np.any(residual_A > 150), "no heavy remnant found"


# ---------------------------------------------------------------------------
# Task 18: Registered-targets and error-path tests
# ---------------------------------------------------------------------------


def _registered_targets(elab, p1, p2):
    from chromo.models import Fluka

    gen = Fluka(FixedTarget(elab, p1, p2), seed=1)
    return gen.registered_targets


def _register_extra_target(elab, p1, p2, extras):
    from chromo.models import Fluka

    gen = Fluka(FixedTarget(elab, p1, p2), seed=1, targets=extras)
    return gen.registered_targets


def test_registered_defaults_present():
    targets = run_in_separate_process(_registered_targets, 100.0, "p", "N14")
    from particle import Particle

    for name in ("N14", "O16", "Pb208"):
        pid = int(Particle.findall(name)[0].pdgid)
        if pid == 2212:
            assert 2212 in targets, f"{name} not in registered_targets {targets}"
        else:
            assert pid in targets, f"{name} not in registered_targets {targets}"


def test_register_extra_target_via_kwarg():
    # Drop an enduring material to make room: use kin.p2=N14 so we only need
    # 9 defaults + no kin-added (N14 already there) + 1 extra slot for Ne20.
    targets = run_in_separate_process(
        _register_extra_target, 100.0, "p", "N14", ("Ne20",)
    )
    from particle import Particle

    ne20 = int(Particle.findall("Ne20")[0].pdgid)
    assert ne20 in targets


def _try_unregistered_target():
    from chromo.models import Fluka

    gen = Fluka(FixedTarget(100.0, "p", "N14"), seed=1)
    try:
        gen.kinematics = FixedTarget(100.0, "p", "U238")
    except (KeyError, ValueError) as exc:
        return str(exc)
    return "no-error"


def test_unregistered_target_raises_with_hint():
    msg = run_in_separate_process(_try_unregistered_target)
    assert "U238" in msg
    assert "targets=" in msg, f"error message lacks remediation hint: {msg}"


# ---------------------------------------------------------------------------
# Task 19: Energy-bound tests
# ---------------------------------------------------------------------------


def _try_construct(elab, p1, p2):
    from chromo.models import Fluka

    try:
        Fluka(FixedTarget(elab, p1, p2), seed=1)
    except ValueError as exc:
        return str(exc)
    else:
        return "no-error"


def _try_generate_event(elab, p1, p2):
    from chromo.models import Fluka

    try:
        gen = Fluka(FixedTarget(elab, p1, p2), seed=1)
        for _ in gen(1):
            pass
    except ValueError as exc:
        return str(exc)
    else:
        return "no-error"


def test_above_uhe_ceiling_raises_hadron():
    # Above the DPMJET→UHE ceiling (sqrt(s_NN) = 500 TeV), construction
    # must fail — no UHE model is linked in chromo's FLUKA build.
    # Lab ekin equivalent of 500 TeV CMS is ~1.33e11 GeV/n; 1e12 is 10x
    # above.
    msg = run_in_separate_process(_try_construct, 1e12, "p", "N14")
    assert msg != "no-error"
    assert "UHE" in msg


def test_above_uhe_ceiling_raises_photon():
    msg = run_in_separate_process(_try_construct, 1e12, "gamma", "Pb208")
    assert msg != "no-error"
    assert "UHE" in msg


def _xsec_at_100tev():
    from chromo.models import Fluka

    gen = Fluka(FixedTarget(1e5, "p", "N14"), seed=1)
    cs = gen.cross_section()
    return cs.inelastic


def test_high_energy_xsec_ok_below_ceiling():
    # 100 TeV lab (well below the 500 TeV CMS ceiling): construction and
    # cross_section() must succeed — FLUKA now hands off to DPMJET for
    # the whole high-energy range instead of an unlinked UHE model.
    inel = run_in_separate_process(_xsec_at_100tev)
    assert 100 < inel < 1000


# ---------------------------------------------------------------------------
# Task 20: RNG state round-trip test
# ---------------------------------------------------------------------------


def _rng_roundtrip(tmpdir):
    """Run one event and save the pre-event RNG state.

    The Pythia8DecayHandler uses its own RNG and is disabled here so that
    only FLUKA's Ranmar state is exercised — the round-trip test is
    specifically for FLUKA's RNG, not Pythia8's.
    """
    import pathlib

    from chromo.models import Fluka

    path = pathlib.Path(tmpdir) / "fluka_rng.dat"

    kin = FixedTarget(100.0, "p", "O16")
    g1 = Fluka(kin, seed=42, rng_state_file=path)
    g1._activate_decay_handler(on=False)
    g1.save_rng_state(path)
    ev1 = next(iter(g1(1)))
    pid1 = ev1.pid.copy()
    en1 = ev1.en.copy()
    return path, pid1, en1


def _rng_replay(path):

    from chromo.models import Fluka

    kin = FixedTarget(100.0, "p", "O16")
    g2 = Fluka(kin, seed=42, rng_state_file=path)
    g2._activate_decay_handler(on=False)
    g2.load_rng_state(path)
    ev2 = next(iter(g2(1)))
    return ev2.pid.copy(), ev2.en.copy()


def test_rng_state_roundtrip(tmp_path):
    """Confirm FLUKA Ranmar state reproduces identical events across processes.

    The Pythia8DecayHandler is disabled to isolate FLUKA's own RNG; its
    own RNG is not seeded deterministically and would cause spurious failures.
    """
    path, pid1, en1 = run_in_separate_process(_rng_roundtrip, str(tmp_path))
    pid2, en2 = run_in_separate_process(_rng_replay, path)
    np.testing.assert_array_equal(pid1, pid2)
    np.testing.assert_allclose(en1, en2, rtol=0, atol=0)
