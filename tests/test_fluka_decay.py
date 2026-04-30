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
