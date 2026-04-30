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
