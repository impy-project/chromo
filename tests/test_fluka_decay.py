"""Tests for chromo's FLUKA radioactive-decay interface."""

import os
import sys

import pytest

from tests.util import run_in_separate_process

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


def _import_smoke():
    from chromo.models import _fluka

    assert hasattr(_fluka, "chromo_dcy_init")
    return True


def test_chromo_dcy_init_symbol_exists():
    assert run_in_separate_process(_import_smoke) is True


def _call_init():
    from chromo.models import _fluka

    _fluka.chromo_dcy_init()
    return True


def test_chromo_dcy_init_runs():
    assert run_in_separate_process(_call_init) is True


def _lookup_u238():
    from chromo.models import _fluka

    _fluka.chromo_dcy_init()
    found, t12, exm, _jsp, _jpt = _fluka.chromo_dcy_lookup(238, 92, 0)
    return int(found), float(t12), float(exm)


def test_dcy_lookup_u238():
    found, t12, exm = run_in_separate_process(_lookup_u238)
    assert found == 1
    assert abs(t12 - 1.41e17) / 1.41e17 < 0.01  # ~4.5 Gyr
    assert abs(exm - 47.31) < 0.05  # MeV


def _lookup_invalid():
    from chromo.models import _fluka

    _fluka.chromo_dcy_init()
    found, _t12, _exm, _jsp, _jpt = _fluka.chromo_dcy_lookup(2, 50, 0)
    return int(found)


def test_dcy_lookup_invalid_returns_zero():
    assert run_in_separate_process(_lookup_invalid) == 0


def _catalog_size():
    import numpy as np

    from chromo.models import _fluka

    _fluka.chromo_dcy_init()

    max_n = 5500
    a = np.zeros(max_n, dtype=np.int32)
    z = np.zeros(max_n, dtype=np.int32)
    m = np.zeros(max_n, dtype=np.int32)
    t12 = np.zeros(max_n, dtype=np.float64)
    exm = np.zeros(max_n, dtype=np.float64)
    jsp = np.zeros(max_n, dtype=np.int32)
    jpt = np.zeros(max_n, dtype=np.int32)

    n = _fluka.chromo_dcy_catalog(max_n, a, z, m, t12, exm, jsp, jpt)
    return int(n), int(a[:n].max()), int(z[:n].max())


def test_dcy_catalog_size():
    n, a_max, z_max = run_in_separate_process(_catalog_size)
    assert 4000 <= n <= 5500
    assert 200 < a_max <= 295
    assert 80 < z_max <= 110


def _channels_cs137():
    import numpy as np

    from chromo.models import _fluka

    _fluka.chromo_dcy_init()
    max_ch = 8
    kind = np.zeros(max_ch, dtype=np.int32)
    br = np.zeros(max_ch, dtype=np.float64)
    da = np.zeros(max_ch, dtype=np.int32)
    dz = np.zeros(max_ch, dtype=np.int32)
    dm = np.zeros(max_ch, dtype=np.int32)
    qv = np.zeros(max_ch, dtype=np.float64)

    n = _fluka.chromo_dcy_channels(137, 55, 0, max_ch, kind, br, da, dz, dm, qv)
    return int(n), kind[:n].tolist(), br[:n].tolist()


def test_dcy_channels_cs137():
    n, kinds, brs = run_in_separate_process(_channels_cs137)
    assert n >= 1
    # All recorded channels should be Beta- (kind=2) for Cs-137
    assert all(k == 2 for k in kinds), kinds
    assert abs(sum(brs) - 1.0) < 0.01


def _channels_cs137_saturated():
    """Caller passes max_ch=1 even though Cs-137 has 2 channels.
    Ensures the BR slot is not corrupted by a write past the cap.
    """
    import numpy as np

    from chromo.models import _fluka

    _fluka.chromo_dcy_init()
    max_ch = 1
    kind = np.zeros(max_ch, dtype=np.int32)
    br = np.zeros(max_ch, dtype=np.float64)
    da = np.zeros(max_ch, dtype=np.int32)
    dz = np.zeros(max_ch, dtype=np.int32)
    dm = np.zeros(max_ch, dtype=np.int32)
    qv = np.zeros(max_ch, dtype=np.float64)
    n = _fluka.chromo_dcy_channels(137, 55, 0, max_ch, kind, br, da, dz, dm, qv)
    return int(n), float(br[0]), int(kind[0]), int(da[0]), int(dz[0])


def test_dcy_channels_saturation_preserves_first_slot():
    n, br0, kind0, da0, dz0 = run_in_separate_process(_channels_cs137_saturated)
    assert n == 1
    # Cs-137 first channel: B- to Ba-137m1, BR ~0.947
    assert kind0 == 2
    assert (da0, dz0) == (137, 56)
    assert abs(br0 - 0.947) < 0.005


def _alloc_line_buffers(max_l):
    import numpy as np

    return (
        np.zeros(max_l, dtype=np.float64),  # br
        np.zeros(max_l, dtype=np.float64),  # e_mev
        np.zeros(max_l, dtype=np.int32),  # nlev
        np.zeros(max_l, dtype=np.int32),  # lpos
    )


def _gamma_lines_ba137m():
    from chromo.models import _fluka

    _fluka.chromo_dcy_init()
    max_l = 64
    br, e_mev, nlev, lpos = _alloc_line_buffers(max_l)
    # The famous 661.66 keV line associated with Cs-137 sources is in
    # FLUKA's table under the metastable daughter Ba-137m (im=1).
    n = _fluka.chromo_dcy_lines(137, 56, 1, 1, max_l, br, e_mev, nlev, lpos)
    return int(n), e_mev[:n].tolist(), br[:n].tolist()


def test_dcy_lines_ba137m_gamma():
    n, energies, brs = run_in_separate_process(_gamma_lines_ba137m)
    assert n >= 1
    # Ba-137m de-excites via the famous 661.66 keV gamma (the line that
    # makes Cs-137/Ba-137m a calibration source).  BR ~0.90 in FLUKA.
    matches = [(e, br) for e, br in zip(energies, brs) if abs(e - 0.66166) < 0.001]
    assert matches, energies
    assert abs(matches[0][1] - 0.899) < 0.01, matches


def _alpha_lines_u238():
    from chromo.models import _fluka

    _fluka.chromo_dcy_init()
    max_l = 64
    br, e_mev, nlev, lpos = _alloc_line_buffers(max_l)
    # U-238 alpha decays. kind=2 (alpha), 3 doubles per record.
    n = _fluka.chromo_dcy_lines(238, 92, 0, 2, max_l, br, e_mev, nlev, lpos)
    return int(n), e_mev[:n].tolist(), br[:n].tolist()


def test_dcy_lines_u238_alpha():
    n, energies, brs = run_in_separate_process(_alpha_lines_u238)
    assert n >= 1
    # Dominant U-238 alpha line is around 4.20 MeV.
    assert any(abs(e - 4.20) < 0.05 for e in energies), energies
    # All BRs are valid probabilities.
    assert all(0.0 <= br <= 1.0 for br in brs), brs


def _beta_spectra_cs137():
    from chromo.models import _fluka

    _fluka.chromo_dcy_init()
    max_l = 64
    br, e_mev, nlev, lpos = _alloc_line_buffers(max_l)
    # Cs-137 has two beta- branches (endpoints ~0.514, 1.176 MeV).
    n = _fluka.chromo_dcy_lines(137, 55, 0, 4, max_l, br, e_mev, nlev, lpos)
    return int(n), e_mev[:n].tolist(), br[:n].tolist(), lpos[:n].tolist()


def test_dcy_lines_cs137_beta():
    n, endpoints, _brs, lpos = run_in_separate_process(_beta_spectra_cs137)
    assert n >= 1
    # Endpoints around 0.514 MeV (94.7%) and 1.176 MeV (5.3%) for Cs-137.
    assert any(abs(e - 0.514) < 0.02 for e in endpoints), endpoints
    assert any(abs(e - 1.176) < 0.02 for e in endpoints), endpoints
    # Cs-137 is pure beta-, no positrons.
    assert all(p == 0 for p in lpos), lpos


def _sample_cs137_one():
    from chromo.models import _fluka

    _fluka.chromo_dcy_init()
    success, kdcy, ilv = _fluka.chromo_dcy_sample(137, 55, 0)
    return bool(success), int(kdcy), int(ilv), int(_fluka.hepevt.nhep)


def test_dcy_sample_cs137_succeeds():
    ok, kdcy, _ilv, nhep = run_in_separate_process(_sample_cs137_one)
    assert ok
    assert kdcy == 2  # B-
    assert nhep >= 2  # at least e- + nubar


def _scaffold_imports():
    from chromo.models.fluka_decay import (  # noqa: F401
        STABLE_DEFAULT,
        DecayChannel,
        DecayLine,
        FlukaDecay,
        FlukaIsotope,
    )

    return True


def test_module_scaffold_imports():
    assert run_in_separate_process(_scaffold_imports) is True


def _dataclass_smoke():
    from chromo.models.fluka_decay import DecayChannel, DecayLine

    ch = DecayChannel(
        name="B-",
        branching=1.0,
        daughter_A=137,
        daughter_Z=56,
        daughter_m=1,
        q_value=0.514,
    )
    ln = DecayLine(energy=0.66166, branching=0.85, end_level=0, is_positron=False)
    return ch.name, ch.daughter_A, ln.energy


def test_dataclass_smoke():
    name, a, e = run_in_separate_process(_dataclass_smoke)
    assert name == "B-"
    assert a == 137
    assert abs(e - 0.66166) < 1e-6


def _make_isotope():
    from chromo.models.fluka_decay import FlukaIsotope

    iso = FlukaIsotope(
        owner=None,
        A=137,
        Z=55,
        m=0,
        t_half=9.49e8,
        mass_excess=-86.546,
        symbol="Cs",
        j_spin=7,
        j_parity=1,
    )
    return (iso.A, iso.Z, iso.m, iso.symbol, iso.short())


def test_isotope_init_and_short():
    A, Z, m, sym, short = run_in_separate_process(_make_isotope)
    assert A == 137 and Z == 55 and m == 0 and sym == "Cs"
    assert "Cs137" in short
    assert "9.49" in short or "9.490" in short


def _construct_and_lookup():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay(seed=42)
    iso = dcy.lookup("Cs137")
    return (iso.A, iso.Z, iso.m, iso.symbol, round(iso.t_half / 1e8, 1))


def test_fluka_decay_construct_and_lookup_cs137():
    A, Z, m, sym, t12_e8 = run_in_separate_process(_construct_and_lookup)
    assert (A, Z, m) == (137, 55, 0)
    assert sym == "Cs"
    assert 9.0 < t12_e8 < 10.0


def _lookup_isomer():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    iso = dcy.lookup("Tc99m")
    return iso.A, iso.Z, iso.m, round(iso.t_half, 0)


def test_fluka_decay_lookup_isomer():
    A, Z, m, t12 = run_in_separate_process(_lookup_isomer)
    assert (A, Z, m) == (99, 43, 1)
    assert 20000 <= t12 <= 23000  # ~ 6.0 hr ≈ 21600 s


def _lookup_tuple_and_missing():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    iso_t = dcy.lookup(238, 92, 0)
    iso_none = dcy.lookup(2, 50, 0)
    return iso_t.A, iso_t.Z, iso_none


def test_fluka_decay_lookup_tuple_and_missing():
    A, Z, none_ = run_in_separate_process(_lookup_tuple_and_missing)
    assert (A, Z) == (238, 92)
    assert none_ is None


def _catalog_full_and_filter():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    full = dcy.catalog()
    long_lived = dcy.catalog(t_half_min=1e10)
    short_only = dcy.catalog(t_half_max=1.0, t_half_min=1e-3)
    actinides = dcy.catalog(z_min=89, z_max=99)
    return len(full), len(long_lived), len(short_only), len(actinides)


def test_catalog_filters():
    n, n_long, n_short, n_act = run_in_separate_process(_catalog_full_and_filter)
    assert 4000 <= n <= 5500
    assert 200 <= n_long <= 600  # plenty of T1/2 > 1e10 s entries
    assert n_short > 0
    assert n_act > 0


def _catalog_contains_known():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    full = dcy.catalog()
    azm = {(i.A, i.Z, i.m) for i in full}
    return (
        (238, 92, 0) in azm,
        (137, 55, 0) in azm,
        (14, 6, 0) in azm,
        (99, 43, 1) in azm,
    )


def test_catalog_contains_known_isotopes():
    found = run_in_separate_process(_catalog_contains_known)
    assert all(found), found


def _cs137_channels():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    iso = dcy.lookup("Cs137")
    chs = iso.channels
    return len(chs), chs[0].name, chs[0].daughter_A, chs[0].daughter_Z


def test_isotope_channels_cs137():
    n, name, dA, dZ = run_in_separate_process(_cs137_channels)
    assert n >= 1
    assert name == "B-"
    assert (dA, dZ) == (137, 56)


def _u238_alpha():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    iso = dcy.lookup("U238")
    return [(c.name, c.branching, c.daughter_A, c.daughter_Z) for c in iso.channels]


def test_isotope_channels_u238_has_alpha():
    chs = run_in_separate_process(_u238_alpha)
    assert any(c[0] == "alpha" and c[2] == 234 and c[3] == 90 for c in chs)


def _cs137_gamma_and_str():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    _iso = dcy.lookup("Ba137")  # Note: 661.66 keV gamma is on Ba-137m
    iso_m = dcy.lookup("Ba137m")
    gammas = iso_m.gamma_lines
    s = str(iso_m)
    return [(round(g.energy, 5), round(g.branching, 4)) for g in gammas], s


def test_gamma_lines_and_str():
    gammas, s = run_in_separate_process(_cs137_gamma_and_str)
    energies = [e for e, _b in gammas]
    assert any(abs(e - 0.66166) < 0.001 for e in energies), gammas
    # str(iso) should contain identifier + gamma energy
    assert "Ba137m" in s or "Ba137" in s
    assert "0.66" in s or "0.6617" in s
