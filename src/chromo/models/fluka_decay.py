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

    @property
    def channels(self) -> tuple[DecayChannel, ...]:
        """Decay channels (lazy fetch on first access)."""
        if self._channels is None:
            self._channels = self._owner._fetch_channels(self.A, self.Z, self.m)
        return self._channels

    def _line_list(self, kind_int: int, kind_str: str) -> tuple[DecayLine, ...]:
        if kind_str in self._lines:
            return self._lines[kind_str]
        result = self._owner._fetch_lines(self.A, self.Z, self.m, kind_int)
        self._lines[kind_str] = result
        return result

    @property
    def gamma_lines(self) -> tuple[DecayLine, ...]:
        return self._line_list(1, "gamma")

    @property
    def alpha_lines(self) -> tuple[DecayLine, ...]:
        return self._line_list(2, "alpha")

    @property
    def ce_lines(self) -> tuple[DecayLine, ...]:
        return self._line_list(3, "ce")

    @property
    def beta_spectra(self) -> tuple[DecayLine, ...]:
        return self._line_list(4, "beta")

    def __str__(self) -> str:
        m_tag = "" if self.m == 0 else f"m{self.m}"
        head = (
            f"Isotope {self.symbol}{self.A}{m_tag}  "
            f"(A={self.A}, Z={self.Z}, m={self.m})\n"
            f"  T1/2     = {self.t_half:.3e} s\n"
            f"  ExMass   = {self.mass_excess:.4f} MeV\n"
            f"  J        = {self.j_spin / 2:.1f}  parity = "
            f"{'+' if self.j_parity > 0 else '-' if self.j_parity < 0 else '?'}"
        )
        rows = []
        for c in self.channels:
            d = (
                f"-> {c.daughter_A:3d}/{c.daughter_Z:3d}/m{c.daughter_m}"
                if c.daughter_A >= 0
                else "-> (no single daughter)"
            )
            rows.append(
                f"    {c.name:<5s}  BR={c.branching * 100:7.3f}%  "
                f"{d}  Q={c.q_value:.4f} MeV"
            )
        chan = "\n  Channels:\n" + "\n".join(rows) if rows else ""

        def _fmt_lines(label, lines, n_max=10):
            if not lines:
                return ""
            out = [f"\n  {label} ({len(lines)}):"]
            out.extend(
                f"    E={ln.energy:9.5f} MeV  BR={ln.branching * 100:7.3f}%"
                for ln in lines[:n_max]
            )
            if len(lines) > n_max:
                out.append(f"    ... ({len(lines) - n_max} more)")
            return "\n".join(out)

        body = (
            chan
            + _fmt_lines("Gamma lines", self.gamma_lines)
            + _fmt_lines("Alpha lines", self.alpha_lines)
            + _fmt_lines("CE/Auger lines", self.ce_lines)
            + _fmt_lines("Beta+/- spectra", self.beta_spectra)
        )
        return head + body


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

        # Eagerly materialise so STABLE_DEFAULT is populated for chain ops.
        self._catalog = self._fetch_full_catalog()
        if not STABLE_DEFAULT:
            _populate_stable_default(self._catalog)

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

    # -- lazy backends used by FlukaIsotope -----------------------------

    def _fetch_channels(self, A: int, Z: int, m: int) -> tuple[DecayChannel, ...]:
        import numpy as np

        max_ch = 8
        kind = np.zeros(max_ch, dtype=np.int32)
        br = np.zeros(max_ch, dtype=np.float64)
        da = np.zeros(max_ch, dtype=np.int32)
        dz = np.zeros(max_ch, dtype=np.int32)
        dm = np.zeros(max_ch, dtype=np.int32)
        qv = np.zeros(max_ch, dtype=np.float64)
        n = self._lib.chromo_dcy_channels(A, Z, m, max_ch, kind, br, da, dz, dm, qv)
        return tuple(
            DecayChannel(
                name=_CHANNEL_NAMES.get(int(kind[i]), "other"),
                branching=float(br[i]),
                daughter_A=int(da[i]),
                daughter_Z=int(dz[i]),
                daughter_m=int(dm[i]),
                q_value=float(qv[i]),
            )
            for i in range(int(n))
        )

    def _fetch_lines(self, A: int, Z: int, m: int, kind: int) -> tuple[DecayLine, ...]:
        import numpy as np

        max_l = 256
        br = np.zeros(max_l, dtype=np.float64)
        e_mev = np.zeros(max_l, dtype=np.float64)
        nlev = np.zeros(max_l, dtype=np.int32)
        lpos = np.zeros(max_l, dtype=np.int32)
        n = self._lib.chromo_dcy_lines(A, Z, m, kind, max_l, br, e_mev, nlev, lpos)
        return tuple(
            DecayLine(
                energy=float(e_mev[i]),
                branching=float(br[i]),
                end_level=int(nlev[i]),
                is_positron=bool(lpos[i]),
            )
            for i in range(int(n))
        )

    # -- inclusive event sampling --------------------------------------

    # FLUKA's RDDCAY tags genuinely stable nuclides with T1/2 = 1e38 s.
    _STABLE_T_HALF = 1e38

    def __call__(self, parent, n: int):
        """Yield ``n`` correlated decay events for ``parent``.

        Parameters
        ----------
        parent : str | tuple | int
            Same forms accepted by ``lookup`` (name, ``(A, Z, m)``, or
            PDG ion code).
        n : int
            Number of events to generate.

        Yields
        ------
        EventData
            One event per yield, with FLUKA's products in
            ``ev.pid / px / py / pz / en / mass``.

        Raises
        ------
        ValueError
            If ``parent`` is unknown to FLUKA's table or stable.
        """
        # Resolve parent first so we validate name and catch stables before
        # the (potentially expensive) sampling loop starts.
        if isinstance(parent, tuple):
            iso = self.lookup(*parent)
        else:
            iso = self.lookup(parent)
        if iso is None:
            msg = f"FLUKA has no decay data for parent {parent!r}."
            raise ValueError(msg)
        if iso.t_half >= self._STABLE_T_HALF:
            msg = (
                f"{parent!r} is stable in FLUKA's table "
                f"(T1/2 = {iso.t_half:.3e} s); cannot decay-sample."
            )
            raise ValueError(msg)

        for _ in range(int(n)):
            yield self._sample_one(iso.A, iso.Z, iso.m)

    def _sample_one(self, A: int, Z: int, m: int):
        """One SPDCEV call -> one EventData built from HEPEVT.

        FlukaDecay has no kinematics/frame, so we bypass FlukaEvent and
        build an EventData directly from the HEPEVT common block populated
        by FLLHEP.  Falls back to an empty EventData on persistent failure
        (logged via ``chromo.util.info``).
        """

        from chromo.util import info as _info

        ok = False
        for _attempt in range(2):
            ok_int, _kdcy, _ilv = self._lib.chromo_dcy_sample(A, Z, m)
            if ok_int:
                ok = True
                break
        if not ok:
            _info(0, f"chromo_dcy_sample failed for ({A},{Z},{m}) twice")

        return self._hepevt_to_event_data()

    def _hepevt_to_event_data(self):
        """Snapshot the HEPEVT common block into a fresh EventData."""
        import numpy as np

        from chromo.common import EventData

        evt = self._lib.hepevt
        nhep = int(evt.nhep)
        sel = slice(0, nhep)

        pid = np.array(evt.idhep[sel], dtype=np.int32, copy=True)
        status = np.array(evt.isthep[sel], dtype=np.int32, copy=True)
        phep = np.asarray(evt.phep[:, sel], dtype=np.float64)
        vhep = np.asarray(evt.vhep[:, sel], dtype=np.float64)
        # phep rows: 0=px, 1=py, 2=pz, 3=E, 4=m. vhep rows: 0=vx, ...
        px = np.array(phep[0], dtype=np.float64, copy=True)
        py = np.array(phep[1], dtype=np.float64, copy=True)
        pz = np.array(phep[2], dtype=np.float64, copy=True)
        en = np.array(phep[3], dtype=np.float64, copy=True)
        mass = np.array(phep[4], dtype=np.float64, copy=True)
        vx = np.array(vhep[0], dtype=np.float64, copy=True)
        vy = np.array(vhep[1], dtype=np.float64, copy=True)
        vz = np.array(vhep[2], dtype=np.float64, copy=True)
        vt = np.array(vhep[3], dtype=np.float64, copy=True)

        # Charge: reuse FLUKA's helper for ordinary particles, recover Z
        # for nuclei from the PDG ion code (10LZZZAAAI) — same logic as
        # FlukaEvent._get_charge but without an MCEvent wrapper.
        charge = self._lib.charge_from_pdg_arr(pid).astype(np.float64)
        absp = np.abs(pid.astype(np.int64))
        is_nucleus = absp >= 1_000_000_000
        z_from_code = ((absp // 10000) % 1000).astype(np.float64)
        charge = np.where(is_nucleus, np.sign(pid) * z_from_code, charge)

        mothers = np.array(evt.jmohep[:, sel].T, dtype=np.int32, copy=True)
        daughters = np.array(evt.jdahep[:, sel].T, dtype=np.int32, copy=True)

        return EventData(
            generator=("FLUKA-decay", self._lib.__name__.split(".")[-1]),
            kin=None,
            nevent=int(evt.nevhep),
            impact_parameter=float("nan"),
            n_wounded=(0, 0),
            production_cross_section=float("nan"),
            pid=pid,
            status=status,
            charge=charge,
            px=px,
            py=py,
            pz=pz,
            en=en,
            m=mass,
            vx=vx,
            vy=vy,
            vz=vz,
            vt=vt,
            mothers=mothers,
            daughters=daughters,
        )

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


# --------------------------------------------------------------------- #
# Stable-default set + DecayChainHandler                                #
# --------------------------------------------------------------------- #

# Lepton/photon/light-nucleon PDG ids that anchor the chain.
_BASE_STABLE_PDG = {
    11,
    -11,  # e-, e+
    12,
    -12,  # nu_e, anu_e
    13,
    -13,  # mu-, mu+
    14,
    -14,  # nu_mu, anu_mu
    16,
    -16,  # nu_tau, anu_tau
    22,  # gamma
    2212,  # proton
    2112,  # neutron
}


def _populate_stable_default(catalog) -> None:
    """Add to STABLE_DEFAULT every (A,Z,m) entry with T1/2 = 1e38."""
    STABLE_DEFAULT.clear()
    STABLE_DEFAULT.update(_BASE_STABLE_PDG)
    for iso in catalog:
        if iso.t_half >= 1e38:
            pdg = 1_000_000_000 + iso.Z * 10_000 + iso.A * 10 + iso.m
            STABLE_DEFAULT.add(pdg)


class DecayChainHandler:
    """Recursive decay-chain post-processor.

    For each event passed to ``expand()``, every product whose PDG id is
    *not* in ``final_state`` *and* corresponds to a FLUKA-decayable
    nuclide is sampled via ``SPDCEV`` and replaced with its products.
    Recursion stops when all products are in ``final_state`` or the
    isotope is stable / has no decay data.
    """

    def __init__(
        self,
        owner: FlukaDecay,
        final_state: set[int] | None = None,
        max_depth: int = 20,
        on_max_depth: str = "raise",
    ):
        self._owner = owner
        self.final_state = (
            set(STABLE_DEFAULT) if final_state is None else set(final_state)
        )
        self.max_depth = int(max_depth)
        if on_max_depth not in {"raise", "warn"}:
            msg = f"on_max_depth must be 'raise' or 'warn', got " f"{on_max_depth!r}"
            raise ValueError(msg)
        self.on_max_depth = on_max_depth

    def expand(self, event):
        """Expand chained decays until all products in final_state."""
        import warnings

        import numpy as np

        owner = self._owner
        # Concatenated product arrays + parent indices.
        pid = list(event.pid.tolist())
        px = list(event.px.tolist())
        py = list(event.py.tolist())
        pz = list(event.pz.tolist())
        en = list(event.en.tolist())
        mass = list(event.m.tolist())
        status = (
            list(event.status.tolist()) if event.status is not None else [1] * len(pid)
        )
        mothers_in = event.mothers if event.mothers is not None else None
        parents = [-1] * len(pid)

        # Active queue: only consider status==1 (final-state) products.
        # Status==4 entries are FLUKA history records (parent + 9999 marker)
        # and must not be re-decayed.
        queue: list[tuple[int, int]] = [
            (i, 0) for i in range(len(pid)) if status[i] == 1
        ]

        while queue:
            i, depth = queue.pop()
            this_pid = pid[i]
            if this_pid in self.final_state:
                continue
            if abs(this_pid) < 1_000_000_000:
                # Standard particle but not in final_state: leave as-is
                # (FLUKA's decay tables don't cover it).
                continue
            A = (abs(this_pid) // 10) % 1000
            Z = (abs(this_pid) // 10_000) % 1000
            m = (abs(this_pid) // 100_000_000) % 10
            iso = owner.lookup(A, Z, m)
            if iso is None or iso.t_half >= 1e38:
                continue
            if depth >= self.max_depth:
                msg = (
                    f"DecayChainHandler hit max_depth={self.max_depth} "
                    f"on ({A},{Z},{m})"
                )
                if self.on_max_depth == "raise":
                    raise RuntimeError(msg)
                warnings.warn(msg, stacklevel=2)
                continue

            ok_int, _kdcy, _ilv = owner._lib.chromo_dcy_sample(A, Z, m)
            if not ok_int:
                continue
            # Sampler already wrote to HEPEVT; just snapshot it.
            sub = owner._hepevt_to_event_data()

            sub_status = (
                sub.status.tolist() if sub.status is not None else [1] * len(sub.pid)
            )
            for j in range(len(sub.pid)):
                # Skip FLUKA history records (status==4): the parent and
                # the 9999 sentinel must not re-enter the chain.
                if int(sub_status[j]) == 4:
                    continue
                pid.append(int(sub.pid[j]))
                px.append(float(sub.px[j]))
                py.append(float(sub.py[j]))
                pz.append(float(sub.pz[j]))
                en.append(float(sub.en[j]))
                mass.append(float(sub.m[j]))
                status.append(int(sub_status[j]))
                parents.append(i)
                if int(sub_status[j]) == 1:
                    queue.append((len(pid) - 1, depth + 1))

            # Mark this product as decayed-away by setting status to 0.
            status[i] = 0

        # Drop placeholders (status==0) before storing.
        keep = [k for k, s in enumerate(status) if s != 0]
        new_pid = np.array([pid[k] for k in keep], dtype=np.int32)
        new_px = np.array([px[k] for k in keep], dtype=np.float64)
        new_py = np.array([py[k] for k in keep], dtype=np.float64)
        new_pz = np.array([pz[k] for k in keep], dtype=np.float64)
        new_en = np.array([en[k] for k in keep], dtype=np.float64)
        new_m = np.array([mass[k] for k in keep], dtype=np.float64)
        new_status = np.array([status[k] for k in keep], dtype=np.int32)
        # Recompute charge from the (possibly extended) PDG list.
        std_charge = owner._lib.charge_from_pdg_arr(new_pid).astype(np.float64)
        absp = np.abs(new_pid.astype(np.int64))
        is_nucleus = absp >= 1_000_000_000
        z_from_code = ((absp // 10000) % 1000).astype(np.float64)
        new_charge = np.where(is_nucleus, np.sign(new_pid) * z_from_code, std_charge)
        # Pad/trim vertex arrays so they have one entry per kept particle.
        n_total = len(pid)

        def _resize(name):
            arr = getattr(event, name, None)
            if arr is None:
                return np.zeros(len(keep), dtype=np.float64)
            old_n = len(arr)
            if n_total > old_n:
                arr = np.concatenate([arr, np.zeros(n_total - old_n, dtype=arr.dtype)])
            return np.array([arr[k] for k in keep], dtype=np.float64)

        new_vx = _resize("vx")
        new_vy = _resize("vy")
        new_vz = _resize("vz")
        new_vt = _resize("vt")
        # Re-map mothers indices to the compacted arrays.
        index_map = {old: new for new, old in enumerate(keep)}
        new_mothers = np.full((len(keep), 2), -1, dtype=np.int32)
        for new_idx, old_idx in enumerate(keep):
            p = parents[old_idx]
            if p >= 0:
                new_mothers[new_idx, 0] = index_map.get(p, -1)
            elif mothers_in is not None and old_idx < len(mothers_in):
                # Preserve any pre-existing mother info from the input event.
                m0 = int(mothers_in[old_idx, 0])
                new_mothers[new_idx, 0] = index_map.get(m0, -1) if m0 >= 0 else -1
                if mothers_in.shape[1] > 1:
                    m1 = int(mothers_in[old_idx, 1])
                    new_mothers[new_idx, 1] = index_map.get(m1, -1) if m1 >= 0 else -1

        from chromo.common import EventData

        return EventData(
            generator=event.generator,
            kin=event.kin,
            nevent=event.nevent,
            impact_parameter=event.impact_parameter,
            n_wounded=event.n_wounded,
            production_cross_section=event.production_cross_section,
            pid=new_pid,
            status=new_status,
            charge=new_charge,
            px=new_px,
            py=new_py,
            pz=new_pz,
            en=new_en,
            m=new_m,
            vx=new_vx,
            vy=new_vy,
            vz=new_vz,
            vt=new_vt,
            mothers=new_mothers,
            daughters=None,
        )
