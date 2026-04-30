"""FLUKA radioactive-decay interface.

Exposes FLUKA 2025.1's nuclear-decay tables and `SPDCEV` event sampler
to chromo users via a separate, kinematics-free `FlukaDecay` class.
See ``docs/superpowers/specs/2026-04-30-fluka-decay-interface-design.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

# --------------------------------------------------------------------- #
# Data containers                                                       #
# --------------------------------------------------------------------- #

# FLUKA channel-kind codes returned by chromo_dcy_channels and (a subset
# of these, 1..6 only) by chromo_dcy_sample.  Keep in sync with the
# Fortran wrappers in chromo_fluka.f.
_CHANNEL_NAMES = {
    1: "alpha",
    2: "B-",
    3: "B+",
    4: "EC",
    5: "IT",
    6: "SF",
    7: "B-N",
    8: "B+P",
    9: "B-2N",
    10: "B-3N",
    11: "B-NA",
    12: "other",
}


@dataclass(frozen=True)
class DecayChannel:
    """One radioactive-decay channel from FLUKA's table.

    Attributes
    ----------
    name : str
        Channel label (``"B-"``, ``"alpha"``, ``"EC"``, ...).
    branching : float
        Branching ratio in [0, 1].
    daughter_A, daughter_Z, daughter_m : int
        Daughter nuclide. ``-1`` if the channel has no single daughter
        (e.g. spontaneous fission).
    q_value : float
        Q value in MeV (atomic, computed via FLUKA's ``QRDDCY``).
    """

    name: str
    branching: float
    daughter_A: int
    daughter_Z: int
    daughter_m: int
    q_value: float


@dataclass(frozen=True)
class DecayLine:
    """One discrete line from FLUKA's gamma/alpha/CE/beta tables.

    Attributes
    ----------
    energy : float
        Line energy in MeV. For B+/-, this is the end-point.
    branching : float
        Branching ratio in [0, 1] of the parent decay.
    end_level : int
        For alpha / B+/-, daughter level index reached. ``0`` for gamma / CE.
    is_positron : bool
        For B+/- only - true if positron, false if electron.
    """

    energy: float
    branching: float
    end_level: int
    is_positron: bool


# Default "stop the chain here" set for `FlukaDecay.chain`. Populated on
# first FlukaDecay instantiation by `_compute_stable_default`.
STABLE_DEFAULT: set[int] = set()


# --------------------------------------------------------------------- #
# FlukaIsotope and FlukaDecay placeholders - filled in later tasks      #
# --------------------------------------------------------------------- #


class FlukaIsotope:
    """Single nuclide record from FLUKA's decay catalogue.

    Filled out across Tasks 8, 11, 12.
    """

    __slots__ = (
        "A",
        "Z",
        "_channels",
        "_lines",
        "_owner",
        "j_parity",
        "j_spin",
        "m",
        "mass_excess",
        "symbol",
        "t_half",
    )


class FlukaDecay:
    """FLUKA radioactive-decay generator and database.

    Filled out across Tasks 9, 10, 13, 14, 15.
    """
