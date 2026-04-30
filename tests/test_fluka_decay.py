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
