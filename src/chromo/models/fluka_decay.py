"""FLUKA radioactive-decay interface.

Exposes FLUKA 2025.1's nuclear-decay tables and `SPDCEV` event sampler
to chromo users via a separate, kinematics-free `FlukaDecay` class.
See ``docs/superpowers/specs/2026-04-30-fluka-decay-interface-design.md``.
"""

from __future__ import annotations

import re
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
# Module-level init guard + helpers                                     #
# --------------------------------------------------------------------- #

# Module-level guard - shared with chromo.models.fluka.Fluka via
# `_mark_fluka_init_done()` (set in fluka.py after STPXYZ in Task 16).
_INITIALIZED = False
_INSTANCE_COUNT = 0
_NAME_RE = re.compile(r"^([A-Za-z]{1,2})(\d{1,3})(?:m(\d*))?$")


def _mark_fluka_init_done():
    """Called by Fluka.__init__ after STPXYZ to share the init state."""
    global _INITIALIZED
    _INITIALIZED = True


def _ensure_init(lib):
    global _INITIALIZED
    if _INITIALIZED:
        return
    lib.chromo_dcy_init()
    _INITIALIZED = True


# Compact Z->symbol table (covers Z=0..110, sufficient for FLUKA's
# experimental decay catalogue):
_Z_FALLBACK = {
    0: "n",
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    21: "Sc",
    22: "Ti",
    23: "V",
    24: "Cr",
    25: "Mn",
    26: "Fe",
    27: "Co",
    28: "Ni",
    29: "Cu",
    30: "Zn",
    31: "Ga",
    32: "Ge",
    33: "As",
    34: "Se",
    35: "Br",
    36: "Kr",
    37: "Rb",
    38: "Sr",
    39: "Y",
    40: "Zr",
    41: "Nb",
    42: "Mo",
    43: "Tc",
    44: "Ru",
    45: "Rh",
    46: "Pd",
    47: "Ag",
    48: "Cd",
    49: "In",
    50: "Sn",
    51: "Sb",
    52: "Te",
    53: "I",
    54: "Xe",
    55: "Cs",
    56: "Ba",
    57: "La",
    58: "Ce",
    59: "Pr",
    60: "Nd",
    61: "Pm",
    62: "Sm",
    63: "Eu",
    64: "Gd",
    65: "Tb",
    66: "Dy",
    67: "Ho",
    68: "Er",
    69: "Tm",
    70: "Yb",
    71: "Lu",
    72: "Hf",
    73: "Ta",
    74: "W",
    75: "Re",
    76: "Os",
    77: "Ir",
    78: "Pt",
    79: "Au",
    80: "Hg",
    81: "Tl",
    82: "Pb",
    83: "Bi",
    84: "Po",
    85: "At",
    86: "Rn",
    87: "Fr",
    88: "Ra",
    89: "Ac",
    90: "Th",
    91: "Pa",
    92: "U",
    93: "Np",
    94: "Pu",
    95: "Am",
    96: "Cm",
    97: "Bk",
    98: "Cf",
    99: "Es",
    100: "Fm",
    101: "Md",
    102: "No",
    103: "Lr",
    104: "Rf",
    105: "Db",
    106: "Sg",
    107: "Bh",
    108: "Hs",
    109: "Mt",
    110: "Ds",
}
_SYM_TO_Z = {v: k for k, v in _Z_FALLBACK.items()}


def _z_to_symbol(z: int) -> str:
    """Element symbol for atomic number Z (or 'n' for free neutron)."""
    return _Z_FALLBACK.get(z, "?")


def _parse_isotope_name(name: str) -> tuple[int, int, int]:
    """Parse 'Cs137', 'Tc99m', 'U238m1' into (A, Z, m).

    'Tc99' -> m=0; 'Tc99m' -> m=1 (1st isomer);
    'U238m2' -> m=2 (2nd isomer).
    """
    match = _NAME_RE.match(name.strip())
    if not match:
        msg = f"Cannot parse isotope name: {name!r}"
        raise ValueError(msg)
    sym, a_str, m_str = match.groups()
    sym = sym[0].upper() + (sym[1:].lower() if len(sym) > 1 else "")
    if sym not in _SYM_TO_Z:
        msg = f"Unknown element symbol: {sym!r}"
        raise ValueError(msg)
    A = int(a_str)
    Z = _SYM_TO_Z[sym]
    if m_str is None:
        m = 0
    elif m_str == "":
        m = 1
    else:
        m = int(m_str)
    return A, Z, m


# --------------------------------------------------------------------- #
# FlukaIsotope and FlukaDecay placeholders - filled in later tasks      #
# --------------------------------------------------------------------- #


class FlukaIsotope:
    """Single nuclide record from FLUKA's decay catalogue.

    Heavy data (decay channels, line lists) is fetched lazily on first
    access via the parent ``FlukaDecay`` instance.  See Tasks 11, 12.
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

    def __init__(
        self, *, owner, A, Z, m, t_half, mass_excess, symbol, j_spin, j_parity
    ):
        self._owner = owner
        self.A = A
        self.Z = Z
        self.m = m
        self.t_half = t_half
        self.mass_excess = mass_excess
        self.symbol = symbol
        self.j_spin = j_spin
        self.j_parity = j_parity
        self._channels = None
        self._lines = {}

    def short(self) -> str:
        """One-line summary."""
        m_tag = "" if self.m == 0 else f"m{self.m}"
        return (
            f"{self.symbol}{self.A}{m_tag} (Z={self.Z}, m={self.m}): "
            f"T1/2={self.t_half:.3e} s, ExM={self.mass_excess:.3f} MeV"
        )

    def __repr__(self):
        return f"<FlukaIsotope {self.short()}>"


class FlukaDecay:
    """FLUKA radioactive-decay generator and database.

    No kinematics are required.  The first instance triggers FLUKA's
    decay-table init (``CMSPPR`` / ``ZEROIN`` / ``NCDTRD`` / ``RDFLUO`` /
    ...).  If a ``Fluka`` instance has already initialised the library
    (via ``STPXYZ``), construction is a no-op.

    Two ``FlukaDecay`` instances in the same process raise
    ``RuntimeError`` (matches the FLUKA singleton).
    """

    def __init__(self, seed: int | None = None):
        global _INSTANCE_COUNT
        if _INSTANCE_COUNT > 0:
            raise RuntimeError(
                "FlukaDecay can only be instantiated once per process "
                "(FLUKA singleton)."
            )
        _INSTANCE_COUNT += 1
        from chromo.models import _fluka

        self._lib = _fluka
        _ensure_init(self._lib)
        if seed is not None:
            self._seed_rng(int(seed))
        self._catalog = None  # filled in Task 10

    @staticmethod
    def _seed_rng(seed: int) -> None:
        """Initialise FLUKA's Ranmar generator with a deterministic seed.

        Mirrors ``Fluka._init_rng``: writes a temporary RNG state file
        via ``init_rng_state(file, lun, seed, ntot, ntot2)``.  ntot=ntot2=0
        means no skipping, logical unit 888 matches the Fluka model.
        """
        import os
        import pathlib
        import tempfile

        from chromo.models import _fluka

        rng_file = (
            pathlib.Path(tempfile.gettempdir())
            / f"fluka_decay_rng_state_{os.getpid()}.dat"
        )
        _fluka.init_rng_state(str(rng_file), 888, int(seed), 0, 0)

    # -- lookup ---------------------------------------------------------

    def lookup(self, *args) -> FlukaIsotope | None:
        """Return the FlukaIsotope for the given (A, Z, m) or name.

        Accepts:

        - ``lookup("Cs137")`` / ``lookup("Tc99m")`` / ``lookup("U238m1")``
        - ``lookup(A, Z, m)``
        - ``lookup(pdg_ion_id)``
        """
        A, Z, m = self._parse_arg(args)
        found, t12, exm, jsp, jpt = self._lib.chromo_dcy_lookup(A, Z, m)
        if not found:
            return None
        return FlukaIsotope(
            owner=self,
            A=A,
            Z=Z,
            m=m,
            t_half=float(t12),
            mass_excess=float(exm),
            symbol=_z_to_symbol(Z) if Z >= 0 else "n",
            j_spin=int(jsp),
            j_parity=int(jpt),
        )

    # -- catalogue ------------------------------------------------------

    def catalog(
        self,
        *,
        t_half_min: float | None = None,
        t_half_max: float | None = None,
        a_min: int | None = None,
        a_max: int | None = None,
        z_min: int | None = None,
        z_max: int | None = None,
    ) -> list[FlukaIsotope]:
        """Return all isotopes/isomers in FLUKA's decay-data table.

        First call materialises the full catalogue (~4 500 entries) by
        calling ``chromo_dcy_catalog``; subsequent calls reuse the cache.

        Parameters
        ----------
        t_half_min, t_half_max : float, optional
            Half-life range in seconds (inclusive).
        a_min, a_max, z_min, z_max : int, optional
            Mass / atomic-number ranges (inclusive).
        """
        if self._catalog is None:
            self._catalog = self._fetch_full_catalog()

        out = self._catalog
        if t_half_min is not None:
            out = [i for i in out if i.t_half >= t_half_min]
        if t_half_max is not None:
            out = [i for i in out if i.t_half <= t_half_max]
        if a_min is not None:
            out = [i for i in out if i.A >= a_min]
        if a_max is not None:
            out = [i for i in out if i.A <= a_max]
        if z_min is not None:
            out = [i for i in out if i.Z >= z_min]
        if z_max is not None:
            out = [i for i in out if i.Z <= z_max]
        return list(out)

    def _fetch_full_catalog(self) -> list[FlukaIsotope]:
        import numpy as np

        max_n = 5500
        a = np.zeros(max_n, dtype=np.int32)
        z = np.zeros(max_n, dtype=np.int32)
        m = np.zeros(max_n, dtype=np.int32)
        t12 = np.zeros(max_n, dtype=np.float64)
        exm = np.zeros(max_n, dtype=np.float64)
        jsp = np.zeros(max_n, dtype=np.int32)
        jpt = np.zeros(max_n, dtype=np.int32)
        n = self._lib.chromo_dcy_catalog(max_n, a, z, m, t12, exm, jsp, jpt)
        return [
            FlukaIsotope(
                owner=self,
                A=int(a[i]),
                Z=int(z[i]),
                m=int(m[i]),
                t_half=float(t12[i]),
                mass_excess=float(exm[i]),
                symbol=_z_to_symbol(int(z[i])),
                j_spin=int(jsp[i]),
                j_parity=int(jpt[i]),
            )
            for i in range(int(n))
        ]

    @staticmethod
    def _parse_arg(args) -> tuple[int, int, int]:
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, str):
                return _parse_isotope_name(arg)
            if isinstance(arg, int):
                # PDG ion code: 10LZZZAAAI
                if arg >= 1_000_000_000:
                    A = (arg // 10) % 1000
                    Z = (arg // 10_000) % 1000
                    m = (arg // 100_000_000) % 10
                    return A, Z, m
                msg = f"Cannot parse single-int arg {arg}"
                raise ValueError(msg)
        if len(args) == 3:
            return int(args[0]), int(args[1]), int(args[2])
        msg = "lookup() expects (name) or (A, Z, m); got " + repr(args)
        raise TypeError(msg)
