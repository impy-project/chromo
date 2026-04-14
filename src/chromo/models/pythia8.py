import math
import warnings
from collections.abc import Collection
from os import environ
from pathlib import Path

import numpy as np
from particle import literals as lp

from chromo.common import CrossSectionData, EventData, MCRun
from chromo.constants import GeV, standard_projectiles
from chromo.kinematics import EventFrame
from chromo.util import Nuclei, _cached_data_dir, name2pdg

# Tabulated average number of inelastic hN collisions per hA collision,
# ported from PythiaCascade.h (Pythia 8.317).  Used by both Cascade and
# Angantyr to compute parametric nuclear cross sections without MC.
_TAB_A = np.array([1, 2, 4, 9, 12, 14, 16, 27, 40, 56, 63, 84, 107, 129, 197, 208])
_TAB_OFFSET = np.array(
    [
        0.0,
        0.0510,
        0.1164,
        0.2036,
        0.2328,
        0.2520,
        0.2624,
        0.3190,
        0.3562,
        0.3898,
        0.3900,
        0.3446,
        0.3496,
        0.3504,
        0.3484,
        0.3415,
    ]
)
_TAB_SLOPE = np.array(
    [
        0.0,
        0.00187,
        0.00496,
        0.0107,
        0.0136,
        0.0152,
        0.0169,
        0.0243,
        0.0314,
        0.0385,
        0.0415,
        0.0506,
        0.0581,
        0.0644,
        0.0806,
        0.0830,
    ]
)
_TAB_SLOPE_LO = np.array(
    [
        0.0,
        0.00361,
        0.00884,
        0.0174,
        0.0210,
        0.0233,
        0.0252,
        0.0340,
        0.0418,
        0.0496,
        0.0524,
        0.0600,
        0.0668,
        0.0727,
        0.0873,
        0.0893,
    ]
)


def _n_coll_avg(A, sigma_hN):
    """Average number of inelastic hN collisions in one hA collision.

    Ported from ``PythiaCascade::nCollAvg`` (PythiaCascade.h, Pythia 8.317).

    Parameters
    ----------
    A : int
        Mass number of the target nucleus (1 ≤ A ≤ 208).
    sigma_hN : float
        Inelastic hadron-nucleon cross section in mb.
    """
    for i, Ai in enumerate(_TAB_A):
        if A == Ai:
            return min(
                1.0 + _TAB_SLOPE_LO[i] * sigma_hN,
                1.0 + _TAB_OFFSET[i] + _TAB_SLOPE[i] * sigma_hN,
            )
        if A < Ai:
            nc1 = min(
                _TAB_SLOPE_LO[i - 1] * sigma_hN,
                _TAB_OFFSET[i - 1] + _TAB_SLOPE[i - 1] * sigma_hN,
            )
            nc2 = min(
                _TAB_SLOPE_LO[i] * sigma_hN,
                _TAB_OFFSET[i] + _TAB_SLOPE[i] * sigma_hN,
            )
            wt1 = (Ai - A) / (Ai - _TAB_A[i - 1])
            return (
                1.0
                + wt1 * (A / _TAB_A[i - 1]) ** (2.0 / 3.0) * nc1
                + (1.0 - wt1) * (A / Ai) ** (2.0 / 3.0) * nc2
            )
    return float("nan")


def _parametric_sigma_hA(A, sigma_hN):
    """Parametric hadron-nucleus inelastic cross section (mb).

    Uses the PythiaCascade formula: sigma(hA) = A * sigma(hN) / nCollAvg(A).
    """
    if A <= 1:
        return sigma_hN
    return A * sigma_hN / _n_coll_avg(A, sigma_hN)


def _merge_cascade_results(results):
    """Concatenate multiple next_coll result tuples, remapping mother/daughter
    indices so they remain valid in the merged particle array.

    Each result is a 13-tuple (pid, status, px, py, pz, en, m, vx, vy, vz, vt,
    mothers, daughters) where mothers/daughters are (N,2) int arrays with
    1-based Pythia indices (0 = no parent/child).
    """
    if len(results) == 1:
        # Caller (_generate) already copies arrays from the C++ buffer.
        return results[0]
    parts = [[] for _ in range(13)]
    offset = 0
    for result in results:
        n = len(result[0])
        for i in range(11):  # scalar arrays: just append
            parts[i].append(result[i])
        for i in (11, 12):  # index arrays: add offset to non-zero entries
            arr = result[i].copy()
            arr[arr > 0] += offset
            parts[i].append(arr)
        offset += n
    return tuple(np.concatenate(parts[i]) for i in range(13))


class PYTHIA8Event(EventData):
    """Wrapper for Pythia8 event stack."""

    def __init__(self, generator):
        pythia = generator._pythia
        event = pythia.event
        super().__init__(
            (generator.name, generator.version),
            generator.kinematics,
            generator.nevents,
            self._get_impact_parameter(pythia),
            self._get_n_wounded(pythia),
            generator._inel_or_prod_cross_section,
            event.pid(),
            event.status(),
            pythia.charge(),
            event.px(),
            event.py(),
            event.pz(),
            event.en(),
            event.m(),
            event.vx(),
            event.vy(),
            event.vz(),
            event.vt(),
            np.maximum(event.mothers() - 1, -1),
            np.maximum(event.daughters() - 1, -1),
        )

    @staticmethod
    def _get_impact_parameter(pythia):
        hi = pythia.info.hiInfo
        if hi is None:
            return np.nan
        return hi.b

    @staticmethod
    def _get_n_wounded(pythia):
        hi = pythia.info.hiInfo
        if hi is None:
            return (0, 0)
        return hi.nPartProj, hi.nPartTarg

    def _prepare_for_hepmc(self):
        model, version = self.generator
        warnings.warn(
            f"{model}-{version}: only part of the history available in HepMC3 event",
            RuntimeWarning,
        )

        # We must apply some workarounds so that HepMC3 conversion and IO works
        # for all models. This should be revisited once the fundamental issues
        # with particle histories have been fixed.

        # to get a valid GenEvent we must
        # 1) select only particles produced after the parton shower
        # 2) connect particles attached to a single beam particle
        # (diffractive events) to the common interaction vertex (1, 2)
        # TODO check if this costs significant amount of time and speed it up if so
        mask = (self.status == 1) | (self.status == 2) | (self.status == 4)
        ev = self[mask]
        mask = (ev.mothers[:, 0] == 0) | (ev.mothers[:, 0] == 1)
        ev.mothers[mask] = (0, 1)
        return ev


class Pythia8(MCRun):
    _name = "Pythia"
    _version = "8.317"
    _library_name = "_pythia8"
    _event_class = PYTHIA8Event
    _frame = EventFrame.CENTER_OF_MASS
    _projectiles = standard_projectiles | {
        p.pdgid
        for p in (
            lp.photon,
            lp.e_plus,
            lp.e_minus,
            # Strange baryons
            lp.Lambda,
            lp.Sigma_plus,
            lp.Sigma_minus,
            lp.Xi_minus,
            lp.Xi_0,
            lp.Omega_minus,
            # Charm mesons
            lp.D_plus,
            lp.D_0,
            lp.D_s_plus,
            # Charm baryons
            lp.Lambda_c_plus,
            # Bottom mesons
            lp.B_plus,
            lp.B_0,
            lp.B_s_0,
            # Bottom baryons
            lp.Lambda_b_0,
        )
    }
    # Nuclear targets/projectiles are NOT supported in the base Pythia8 class.
    # Use Pythia8Angantyr for heavy-ion collisions or Pythia8Cascade for
    # single-collision hadron-nucleus mode.
    _targets = _projectiles | {lp.e_plus.pdgid, lp.e_minus.pdgid}
    _restartable = True
    _data_url = (
        "https://github.com/impy-project/chromo"
        "/releases/download/zipped_data_v1.0/Pythia8_v007.zip"
    )

    def __init__(self, evt_kin, *, seed=None, config=None, banner=True):
        """

        Parameters
        ----------
        config: str or list of str, optional
            Pass standard Pythia configuration here. You can use this to change
            the physics processes which are enabled in Pythia. You can pass a
            single string that is read from a configuration file or a list of
            strings, where each string is a single configuration command.
            If config is not set, 'SoftQCD:inelastic = on' is used to get the
            equivalent of other generators in chromo.
        """

        super().__init__(seed)

        datdir = _cached_data_dir(self._data_url) + "xmldoc"

        # Must delete PYTHIA8DATA from environ if it exists, since it overrides
        # our argument here. When you install Pythia8 with conda, it sets
        # PYTHIA8DATA. If that version does not match this version, readString()
        # or init() may fail.
        if "PYTHIA8DATA" in environ:
            del environ["PYTHIA8DATA"]
        self._pythia = self._lib.Pythia(datdir, banner)

        if config is None:
            if evt_kin.p1 == lp.photon.pdgid and evt_kin.p2 == lp.photon.pdgid:
                self._config = ["PhotonCollision:all = on"]
            else:
                # includes gamma p processes
                self._config = ["SoftQCD:inelastic = on"]
        else:
            self._config = self._parse_config(config)

        # Common settings
        self._config += [
            # use our random seed
            "Random:setSeed = on",
            # Pythia's RANMAR PRNG accepts only seeds smaller than 900_000_000,
            # this may change in the future if they switch to a different PRNG
            f"Random:seed = {self.seed % 900_000_000}",
        ]

        # Add "Print:quiet = on" if no "Print:quiet" is in config
        if not any("Print:quiet" in s for s in self._config):
            self._config.append("Print:quiet = on")

        # must come last
        if evt_kin is None:
            # Decay mode
            self._init_pythia(self._config)
        else:
            self.kinematics = evt_kin
        self._set_final_state_particles()

    def _cross_section(self, kin=None, max_info=False):
        st = self._pythia.info.sigmaTot

        if (self.kinematics.p1.is_lepton) and (self.kinematics.p2.is_lepton):
            return CrossSectionData(
                total=st.sigmaTot,
                inelastic=st.sigmaTot - st.sigmaEl,
                elastic=st.sigmaEl,
            )
        return CrossSectionData(
            total=st.sigmaTot,
            inelastic=st.sigmaTot - st.sigmaEl,
            elastic=st.sigmaEl,
            diffractive_ax=st.sigmaAX,
            diffractive_xb=st.sigmaXB,
            diffractive_xx=st.sigmaXX,
            diffractive_axb=st.sigmaAXB,
        )

    def _set_kinematics(self, kin):
        config = self._config[:]

        # TODO use numpy PRNG instead of Pythia's
        config += [
            # use center-of-mass frame
            "Beams:frameType = 1",
            # do not print progress
            "Next:numberCount = 0",
        ]

        config += [
            f"Beams:idA = {int(kin.p1)}",
            f"Beams:idB = {int(kin.p2)}",
            f"Beams:eCM = {kin.ecm}",
        ]

        self._init_pythia(config)

    def _init_pythia(self, config):
        pythia = self._pythia
        pythia.settings.resetAll()

        for line in config:
            if not pythia.readString(line):
                msg = f"readString({line!r}) failed"
                raise RuntimeError(msg)

        self._final_config = config

        # calling init several times is allowed
        if not pythia.init():
            raise RuntimeError("Pythia8 initialization failed")

    def _set_stable(self, pdgid, stable):
        self._pythia.particleData.mayDecay(pdgid, not stable)

    def _get_stable(self):
        r = set()
        for p in self._pythia.particleData.all():
            # the p.tau0 > 1e-5 cut should not be necessary, but without it,
            # Pythia8 reports nuclei and some BSM particles like tau' as stable
            if p.tau0 > 1e-5 and not p.mayDecay:
                r.add(p.id)
                if p.hasAnti:
                    r.add(-p.id)
        return r

    def _generate(self):
        return self._pythia.next()

    @staticmethod
    def _parse_config(config) -> list[str]:
        # convert config to lines and filter out lines that
        # we override with our settings in _set_kinematics
        result: list[str] = []
        if isinstance(config, str):
            for line in config.split("\n"):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                result.append(line)
        elif isinstance(config, Collection):
            result.extend([item.strip() for item in config])

        ignored = ("Random:", "Beams:", "Next:")

        result2: list[str] = []
        for line in result:
            line = line.strip()
            for ig in ignored:
                if line.startswith(ig):
                    warnings.warn(f"configuration ignored: {line!r}", RuntimeWarning)
                    break
            else:
                result2.append(line)
        return result2

    @property
    def random_state(self):
        """Get Pythia8's random number generator state."""
        return self._pythia.getRndmState()

    @random_state.setter
    def random_state(self, rng_state):
        """Restore Pythia8's random number generator state."""
        self._pythia.setRndmState(rng_state)


class PYTHIA8CascadeEvent(EventData):
    """Event from a single hadron+nucleus inelastic collision via PythiaCascade."""

    def __init__(self, generator):
        cascade = generator._cascade
        pid, status, px, py, pz, en, m, vx, vy, vz, vt, mothers, daughters = (
            generator._last_result
        )
        super().__init__(
            (generator.name, generator.version),
            generator.kinematics,
            generator.nevents,
            np.nan,
            generator._n_wounded,
            generator._inel_or_prod_cross_section,
            pid,
            status,
            cascade.charge(generator._last_result),
            px,
            py,
            pz,
            en,
            m,
            vx,
            vy,
            vz,
            vt,
            np.maximum(mothers - 1, -1),
            np.maximum(daughters - 1, -1),
        )


class Pythia8Cascade(MCRun):
    """Pythia8 single-collision mode via PythiaCascade.

    Simulates inelastic hadron+nucleus collisions using Pythia8's
    PythiaCascade plugin.  Each call to the generator yields one
    h+A inelastic collision event.

    Nuclear projectiles (A > 1) are decomposed into Z protons and
    (A-Z) neutrons each carrying 1/A of the total lab momentum.
    Each nucleon is collided independently; the last successful
    sub-collision is returned as the event.

    Parameters
    ----------
    evt_kin : kinematics object
        Fixed-target kinematics.  p2 specifies the target nucleus.
    seed : int or None
        Random seed (not forwarded to Pythia; PythiaCascade uses its
        own internal RNG seeded from system entropy).
    eKinMin : float
        Minimum kinetic energy in GeV for a projectile to interact.
        Default 0.3 GeV.
    enhanceSDtarget : float
        Enhancement factor for target-side single diffraction in
        sub-collisions after the first.  Default 0.5.
    init_file : str
        Path to the MPI initialisation .cmnd file.  Leave empty to
        use the bundled ``setups/InitDefaultMPI.cmnd``.
    """

    _name = "Pythia8Cascade"
    _version = "8.317"
    _library_name = "_pythia8"
    _event_class = PYTHIA8CascadeEvent
    _frame = EventFrame.FIXED_TARGET
    # PythiaCascade supports any hadron that Pythia's getSigmaTotal can handle.
    # This includes standard projectiles, strange/charm/bottom hadrons, and nuclei.
    # Nuclear projectiles are decomposed into Z protons and (A-Z) neutrons,
    # each carrying 1/A of the total lab momentum.
    _projectiles = (
        standard_projectiles
        | {
            p.pdgid
            for p in (
                # Strange baryons
                lp.Lambda,
                lp.Sigma_plus,
                lp.Sigma_minus,
                lp.Xi_minus,
                lp.Xi_0,
                lp.Omega_minus,
                # Charm mesons
                lp.D_plus,
                lp.D_0,
                lp.D_s_plus,
                # Charm baryons
                lp.Lambda_c_plus,
                # Bottom mesons
                lp.B_plus,
                lp.B_0,
                lp.B_s_0,
                # Bottom baryons
                lp.Lambda_b_0,
            )
        }
        | Nuclei(a_min=2, a_max=208)
    )
    _targets = {name2pdg(x) for x in ("He4", "N14", "O16", "Ar40", "Fe56", "Pb208")}
    _restartable = True
    _data_url = Pythia8._data_url

    def __init__(
        self,
        evt_kin,
        *,
        seed=None,
        eKinMin=0.3,
        enhanceSDtarget=0.5,
        init_file="",
        banner=True,
    ):
        super().__init__(seed)

        datdir = _cached_data_dir(self._data_url) + "xmldoc"
        # PythiaCascade constructs its own internal Pythia objects, which
        # locate XML data via PYTHIA8DATA.  Set it explicitly here.
        environ["PYTHIA8DATA"] = datdir
        if not init_file:
            bundled = Path(datdir).parent / "setups" / "InitDefaultMPI.cmnd"
            src_tree = (
                Path(__file__).parent.parent.parent
                / "cpp/pythia83/share/Pythia8/setups/InitDefaultMPI.cmnd"
            )
            if bundled.exists():
                init_file = str(bundled)
            elif src_tree.exists():
                init_file = str(src_tree)
            else:
                raise FileNotFoundError(
                    "Could not find InitDefaultMPI.cmnd. "
                    "Pass init_file= explicitly or add setups/ to the data bundle."
                )

        self._cascade = self._lib.PythiaCascadeForChromo(banner)
        # rapidDecays=True with smallTau0=1000 gives the default Pythia8
        # decay setup (see PythiaCascade.h line 113-114).
        # slowDecays=True decays mu/pi/K/K0L, matching cosmic-ray
        # conventions where all particles decay at the interaction point.
        if not self._cascade.init(
            eKinMin,
            enhanceSDtarget,
            init_file,
            rapidDecays=True,
            smallTau0=1000.0,
            slowDecays=True,
            listFinalOnly=True,
        ):
            raise RuntimeError(
                "Pythia8Cascade (PythiaCascadeForChromo) initialization failed"
            )

        self._last_result = None
        self._n_wounded = (0, 0)
        self.kinematics = evt_kin
        self._set_final_state_particles()

    def _set_kinematics(self, kin):
        p2 = kin.p2
        self._Z_targ = p2.Z or 1
        self._A_targ = p2.A or 1

        p1 = kin.p1
        A1 = p1.A or 1
        Z1 = p1.Z or 1
        if A1 > 1:
            self._nucleon_ids = [2212] * Z1 + [2112] * (A1 - Z1)
        else:
            self._nucleon_ids = [int(p1)]

        # Per-nucleon lab momentum and energy
        self._plab_per_nuc = kin.plab / A1
        self._elab_per_nuc = kin.elab / A1

    def _generate(self):
        pd = self._cascade.particle_data()
        results = []
        total_n_coll = 0
        for nid in self._nucleon_ids:
            m = pd.findParticle(nid).m0
            pz = self._plab_per_nuc
            e = math.sqrt(pz * pz + m * m)
            result = self._cascade.next_coll(
                nid, 0.0, 0.0, pz, e, m, self._Z_targ, self._A_targ
            )
            if result is not None:
                # Copy arrays now; the next next_coll() overwrites the C++ buffer.
                results.append(tuple(np.array(a, copy=True) for a in result))
                total_n_coll += self._cascade.n_collisions()
        if not results:
            return False
        self._last_result = _merge_cascade_results(results)
        # Convention differs from Glauber models: first element is the total
        # number of inelastic sub-collisions inside the nucleus (from
        # PythiaCascade::nCollisions), second is the number of projectile
        # nucleons that interacted successfully.
        self._n_wounded = (total_n_coll, len(results))
        return True

    def _cross_section(self, kin=None, max_info=False):
        if kin is None:
            kin = self.kinematics
        pd = self._cascade.particle_data()
        nid = self._nucleon_ids[0]
        m = pd.findParticle(nid).m0
        pz = self._plab_per_nuc
        e = math.sqrt(pz * pz + m * m)
        sigma = self._cascade.sigma_hA(nid, 0.0, 0.0, pz, e, m, self._A_targ)
        return CrossSectionData(inelastic=sigma, prod=sigma)

    def _set_stable(self, pdgid, stable):
        self._cascade.set_may_decay(pdgid, not stable)

    @property
    def random_state(self):
        """Get RNG state for both internal Pythia instances."""
        return {"cascade": self._cascade.getRndmState()}

    @random_state.setter
    def random_state(self, rng_state):
        """Restore RNG state for both internal Pythia instances."""
        self._cascade.setRndmState(rng_state["cascade"])


class Pythia8Angantyr(MCRun):
    """Pythia8 Angantyr heavy-ion model for nuclear interactions.

    Uses precomputed Angantyr initialization grids for fast, proper nuclear
    collision modeling with Glauber geometry, impact parameters, and wounded
    nucleon tracking.  The bundled tables cover CMS energies from 20 GeV to
    20 PeV, which covers cosmic-ray lab energies up to ~200 EeV.

    Parameters
    ----------
    evt_kin : kinematics object
        Center-of-mass or fixed-target kinematics.
    seed : int or None
        Random seed.
    config : str or list of str, optional
        Extra Pythia8 configuration strings applied on top of the Angantyr
        defaults.
    banner : bool
        Whether to print the Pythia8 banner.
    """

    _name = "Pythia8Angantyr"
    _version = "8.317"
    _library_name = "_pythia8"
    _event_class = PYTHIA8Event
    _frame = EventFrame.CENTER_OF_MASS
    _ecm_min = 20 * GeV  # precomputed Angantyr tables start at 20 GeV CMS
    _projectiles = standard_projectiles | Nuclei(a_min=2, a_max=208)
    _targets = {
        name2pdg(x)
        for x in (
            "He4",
            "Li6",
            "C12",
            "N14",
            "O16",
            "Ar40",
            "Cu63",
            "Xe129",
            "Au197",
            "Pb208",
        )
    }
    _restartable = True
    _data_url = Pythia8._data_url

    def __init__(self, evt_kin, *, seed=None, config=None, banner=True):
        super().__init__(seed)

        datdir = _cached_data_dir(self._data_url) + "xmldoc"

        if "PYTHIA8DATA" in environ:
            del environ["PYTHIA8DATA"]
        self._pythia = self._lib.Pythia(datdir, banner)

        # Locate Angantyr setup files (prefer bundled in data zip, fall back to source tree)
        setups_candidates = [
            Path(datdir).parent / "setups",
            Path(__file__).parent.parent.parent / "cpp/pythia83/share/Pythia8/setups",
        ]
        setups_dir = None
        for candidate in setups_candidates:
            if (candidate / "AngantyrCascade.cmnd").exists():
                setups_dir = candidate
                break
        if setups_dir is None:
            raise FileNotFoundError(
                "Could not find Angantyr setup files (AngantyrCascade.cmnd). "
                "Ensure the Pythia8 submodule is checked out or add setups/ to the data bundle."
            )
        self._setups_dir = setups_dir
        self._angantyr_cmnd = setups_dir / "AngantyrCascade.cmnd"
        self._init_angantyr_cmnd = setups_dir / "InitDefaultAngantyr.cmnd"

        self._extra_config = [] if config is None else Pythia8._parse_config(config)
        self._banner = banner
        self._initialized = False
        self._internal_cs_call = False

        self.kinematics = evt_kin
        self._set_final_state_particles()

    @staticmethod
    def _load_cmnd_file(pythia, path):
        """Load a .cmnd file line-by-line, skipping comments and include directives."""
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(("!", "include")):
                    continue
                if not pythia.readString(line):
                    msg = f"readString({line!r}) failed"
                    raise RuntimeError(msg)

    _n_glauber_trials = 10000

    def _init_pythia(self, kin, glauber_only=False):
        """Load setup files and apply settings, then call init()."""
        pythia = self._pythia
        pythia.settings.resetAll()

        # Load the three setup files manually in order, skipping `include =`
        # directives since include paths are resolved relative to the data
        # bundle's xmldoc/../ directory, and setups/ is not in the bundle.
        #
        # Order: MPI init → Angantyr init (includes MPI) → Cascade settings
        self._load_cmnd_file(pythia, self._setups_dir / "InitDefaultMPI.cmnd")
        self._load_cmnd_file(pythia, self._init_angantyr_cmnd)
        self._load_cmnd_file(pythia, self._angantyr_cmnd)

        # Common settings; turn off cascadeMode (which AngantyrCascade.cmnd
        # enables for hadronic air-shower cascades) since here we generate
        # standalone events and expect next() to always return an interaction.
        for line in [
            "Random:setSeed = on",
            f"Random:seed = {self.seed % 900_000_000}",
            "Print:quiet = on",
            "Next:numberCount = 0",
            "Beams:frameType = 1",
            "Beams:allowVariableEnergy = on",
            "Angantyr:cascadeMode = off",
        ]:
            if not pythia.readString(line):
                msg = f"readString({line!r}) failed"
                raise RuntimeError(msg)

        if glauber_only:
            pythia.readString("Angantyr:GlauberOnly = on")

        for line in self._extra_config:
            if not pythia.readString(line):
                msg = f"readString({line!r}) failed"
                raise RuntimeError(msg)

        pythia.readString(f"Beams:idA = {int(kin.p1)}")
        pythia.readString(f"Beams:idB = {int(kin.p2)}")
        pythia.readString(f"Beams:eCM = {kin.ecm}")

        if not pythia.init():
            raise RuntimeError("Pythia8Angantyr initialization failed")

    def _has_precomputed_tables(self):
        """Check if the InitDefaultAngantyr.cmnd contains precomputed tables.

        When Init:reuseHeavyIonSigFit is present, Angantyr init is fast and
        we can afford a GlauberOnly pre-init to compute cross sections.
        Without it, init involves expensive fitting and must not be repeated.
        """
        with open(self._init_angantyr_cmnd) as f:
            for line in f:
                line = line.strip()
                if line.startswith(("!", "#")):
                    continue
                if "Init:reuseHeavyIonSigFit" in line:
                    return True
        return False

    def _run_glauber_trials(self, n_trials=None):
        """Run GlauberOnly trials and return the resulting cross sections."""
        if n_trials is None:
            n_trials = self._n_glauber_trials
        hi = self._pythia.info.hiInfo
        hi.glauberReset()
        for _ in range(n_trials):
            self._pythia.next()
        return CrossSectionData(
            total=hi.glauberTot(),
            inelastic=hi.glauberINEL(),
            elastic=hi.glauberEL(),
            prod=hi.glauberINEL(),
        )

    def _compute_parametric_cs(self, kin):
        """Compute cross section using the PythiaCascade parametric formula.

        Uses ``Pythia::getSigmaTotal/Partial`` for the hadron-nucleon cross
        section and the PythiaCascade nCollAvg table for the nuclear
        scaling.  Fast (no event generation) and deterministic.
        """
        p1 = kin.p1
        # For nuclear projectiles, use proton as representative nucleon
        proj_id = 2212 if (p1.is_nucleus and (p1.A or 1) > 1) else int(p1)
        # Pythia::getSigmaTotal/Partial work correctly in Angantyr mode
        # (unlike info.sigmaTot which holds frozen placeholder values).
        sigma_tot = self._pythia.getSigmaTotal(proj_id, 2212, kin.ecm)
        sigma_el = self._pythia.getSigmaPartial(proj_id, 2212, kin.ecm, 2)
        sigma_hN_inel = sigma_tot - sigma_el

        A_targ = kin.p2.A or 1
        sigma = _parametric_sigma_hA(A_targ, sigma_hN_inel)
        return CrossSectionData(inelastic=sigma, prod=sigma)

    def _set_kinematics(self, kin):
        if not self._initialized:
            if self._has_precomputed_tables():
                # Fast path: precomputed tables make init cheap, so we can
                # afford a GlauberOnly init to pre-compute cross sections
                # followed by a real init for event generation.
                self._init_pythia(kin, glauber_only=True)
                self._glauber_cs = self._run_glauber_trials()
                self._init_pythia(kin, glauber_only=False)
            else:
                # Slow path: no precomputed tables, init is expensive.
                self._glauber_cs = None
                self._init_pythia(kin, glauber_only=False)
            self._initialized = True
            self._idA = int(kin.p1)
            self._idB = int(kin.p2)
        else:
            # Use live switching via setBeamIDs — avoids expensive
            # full re-initialization (~5s per switch).
            idA, idB = int(kin.p1), int(kin.p2)
            if idA != self._idA or idB != self._idB:
                if not self._pythia.setBeamIDs(idA, idB):
                    msg = f"setBeamIDs({idA}, {idB}) failed"
                    raise RuntimeError(msg)
            if not self._pythia.setKinematics(kin.ecm):
                msg = f"setKinematics({kin.ecm}) failed"
                raise RuntimeError(msg)
            self._glauber_cs = None
            self._idA = idA
            self._idB = idB
        # Always available immediately — no event generation needed.
        self._parametric_cs = self._compute_parametric_cs(kin)
        # Flag suppresses the user-facing warning during internal
        # kinematics setter calls; cleared after cross_section() runs.
        self._internal_cs_call = True

    def _generate(self):
        return self._pythia.next()

    def _cross_section(self, kin=None, max_info=False):
        if self._glauber_cs is not None:
            return self._glauber_cs
        return self._parametric_cs

    def cross_section(self, kin=None, max_info=False):
        """Cross sections for the current beam/energy configuration.

        By default returns a fast parametric estimate (PythiaCascade formula).
        For full Angantyr Glauber cross sections, call
        :meth:`glauber_cross_section` instead.
        """
        internal = self._internal_cs_call
        self._internal_cs_call = False
        result = super().cross_section(kin, max_info=max_info)
        if not internal and self._glauber_cs is None:
            warnings.warn(
                "Returning parametric (PythiaCascade) cross section. "
                "For Angantyr Glauber cross sections, call "
                "glauber_cross_section(n_trials).",
                RuntimeWarning,
                stacklevel=2,
            )
        return result

    def glauber_cross_section(self, kin=None, n_trials=10000):
        """Compute Angantyr Glauber cross sections by running MC trials.

        Parameters
        ----------
        kin : EventKinematics, optional
            If provided, compute for these kinematics instead of the
            current setup (same convention as ``cross_section(kin)``).
        n_trials : int
            Number of GlauberOnly events to run (default 10 000).

        Returns
        -------
        CrossSectionData
            Total, inelastic, elastic, and production cross sections.
        """
        if not self._has_precomputed_tables():
            raise RuntimeError(
                "glauber_cross_section() requires precomputed Angantyr tables. "
                "Use the default parametric cross_section() instead."
            )
        if kin is not None:
            self._check_kinematics(kin)
        kin = kin or self.kinematics
        restore = self.kinematics
        self._init_pythia(kin, glauber_only=True)
        result = self._run_glauber_trials(n_trials)
        # Re-init for event generation and restore original kinematics
        self._init_pythia(restore, glauber_only=False)
        self._idA = int(restore.p1)
        self._idB = int(restore.p2)
        self._parametric_cs = self._compute_parametric_cs(restore)
        self._glauber_cs = result if kin is restore else None
        return result

    def _set_stable(self, pdgid, stable):
        self._pythia.particleData.mayDecay(pdgid, not stable)

    def _get_stable(self):
        r = set()
        for p in self._pythia.particleData.all():
            if p.tau0 > 1e-5 and not p.mayDecay:
                r.add(p.id)
                if p.hasAnti:
                    r.add(-p.id)
        return r

    @property
    def random_state(self):
        """Get Pythia8's random number generator state."""
        return self._pythia.getRndmState()

    @random_state.setter
    def random_state(self, rng_state):
        """Restore Pythia8's random number generator state."""
        self._pythia.setRndmState(rng_state)
