import warnings
from collections.abc import Collection
from os import environ

import numpy as np
from particle import literals as lp

from chromo.common import CrossSectionData, EventData, MCRun
from chromo.constants import standard_projectiles
from chromo.kinematics import EventFrame
from chromo.util import _cached_data_dir, name2pdg


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
    _version = "8.315"
    _library_name = "_pythia8"
    _event_class = PYTHIA8Event
    _frame = EventFrame.CENTER_OF_MASS
    _projectiles = standard_projectiles | {
        lp.photon.pdgid,
        lp.e_plus.pdgid,
        lp.e_minus.pdgid,
    }
    # Nuclei are supported in principle, but generation is very slow.
    # Support for more nuclei can be added with ParticleData.addParticle.
    _targets = (
        _projectiles
        | {
            name2pdg(x)
            for x in ("He4", "Li6", "C12", "O16", "Cu63", "Xe129", "Au197", "Pb208")
        }
        | {lp.e_plus.pdgid, lp.e_minus.pdgid}
    )
    _restartable = True
    _data_url = (
        "https://github.com/impy-project/chromo"
        "/releases/download/zipped_data_v1.0/Pythia8_v005.zip"
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

        if (kin.p1.A or 0) > 1 or (kin.p2.A or 1) > 1:
            import warnings

            warnings.warn(
                "Support for nuclei is experimental; event generation takes a long time",
                RuntimeWarning,
            )

            # speed-up initialization by not fitting nuclear cross-section
            config.append("HeavyIon:SigFitNGen = 0")
            config.append("HeavyIon:SigFitDefPar = 10.79,1.75,0.30,0.0,0.0,0.0,0.0,0.0")

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
            result2.append(line)
        return result2


class PYTHIA8CascadeEvent(EventData):
    """Wrapper for Pythia8 cascade event stack."""

    def __init__(self, generator):
        """
        Parameters
        ----------
        generator : Pythia8Cascade
            The generator instance
        """
        # Get the current event from the generator
        event = generator.current_event

        super().__init__(
            (generator.name, generator.version),
            generator.kinematics,
            generator.nevents,
            np.nan,  # impact parameter not available in cascade mode
            (0, 0),  # n_wounded not available in cascade mode
            generator._inel_or_prod_cross_section,
            event.pid(),
            event.status(),
            generator._cascade.charge(event),
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

class Pythia8Cascade(Pythia8):
    """
    Pythia8 wrapper using PythiaCascade for hadron-nucleus interactions.

    This class extends the base Pythia8 class to use the PythiaCascade
    functionality for simulating hadron-nucleus collisions at various energies.
    It's particularly useful for cosmic ray simulations and high-energy
    hadron-nucleus interactions.
    """

    _name = "Pythia8Cascade"
    _version = "8.315"
    _library_name = "_pythia8"
    _event_class = PYTHIA8CascadeEvent
    _frame = EventFrame.FIXED_TARGET

    # Support more nucleus types for cascade interactions
    _targets = standard_projectiles | {
        name2pdg(x)
        for x in (
            "He4",
            "Li6",
            "Be9",
            "C12",
            "N14",
            "O16",
            "Al27",
            "Ar40",
            "Fe56",
            "Cu63",
            "Kr84",
            "Ag107",
            "Xe129",
            "Au197",
            "Pb208",
        )
    }

    def __init__(
        self,
        evt_kin,
        *,
        seed=None,
        config=None,
        banner=True,
        max_energy=1e9,
        use_mpi_tables=True,
        mpi_file="pythiaCascade.mpi",
    ):
        """
        Parameters
        ----------
        max_energy : float, optional
            Maximum energy for cascade initialization (GeV). Default: 1e9
        use_mpi_tables : bool, optional
            Whether to use pre-computed MPI tables for acceleration. Default: True
        mpi_file : str, optional
            File name for MPI initialization data. Default: "pythiaCascade.mpi"
        """

        if evt_kin is None:
            raise ValueError("evt_kin is required for Pythia8Cascade")

        # Ensure we're working in fixed target frame for equivalence with 484 example
        if evt_kin.frame != EventFrame.FIXED_TARGET:
            raise ValueError("Pythia8Cascade requires fixed target frame kinematics")

        # Store cascade-specific parameters
        self._max_energy = max_energy
        self._use_mpi_tables = use_mpi_tables
        self._mpi_file = mpi_file

        # Set default config for cascade mode
        if config is None:
            config = ["SoftQCD:inelastic = on"]

        # Initialize base class with kinematics
        super().__init__(evt_kin, seed=seed, config=config, banner=banner) # TODO: No need to initialize Pythia at all.

        # Initialize cascade after base initialization
        self._init_cascade()

        # Initialize current event storage
        self._current_event = None

    def _init_cascade(self):
        """Initialize the PythiaCascade object."""
        from chromo.util import _cached_data_dir

        # Set up data directory for Pythia XML files
        datdir = _cached_data_dir(self._data_url) + "xmldoc"
        environ["PYTHIA8DATA"] = datdir

        # Initialize cascade with appropriate settings
        reuse_mpi = 3 if self._use_mpi_tables else 0
        init_file = self._mpi_file if self._use_mpi_tables else ""

        # Use main424.cmnd if available for acceleration
        try:
            # Try to find main424.cmnd file
            import os

            possible_paths = [
                "/home/anatoli/devel/chromo/src/cpp/pythia83/examples/main424.cmnd",
                "main424.cmnd",
            ]
            for path in possible_paths:
                if os.path.exists(path):
                    reuse_mpi = -1  # Use cmnd file
                    init_file = path
                    break
        except Exception:
            pass

        self._cascade = self._lib.PythiaCascade()

        # Initialize with appropriate parameters
        if not self._cascade.init(
            self._max_energy,
            False,  # listFinal
            True,  # rapidDecays
            0.0,  # smallTau0
            reuse_mpi,
            init_file,
        ):
            raise RuntimeError("PythiaCascade initialization failed")

    def _set_kinematics(self, kin):
        """Set up kinematics for Pythia 8 Cascade mode collision."""
        # Ensure we're working in fixed target frame for equivalence with 484 example
        if kin.frame != EventFrame.FIXED_TARGET:
            raise ValueError("Pythia8Cascade requires fixed target frame kinematics")

    def _cross_section(self, kin=None, max_info=False):
        """Calculate cross sections for hadron-nucleus interaction."""
        if kin is None:
            kin = self.kinematics

        # Check if cascade is initialized
        if not hasattr(self, "_cascade") or self._cascade is None:
            warnings.warn(
                "Cascade not initialized, returning zero cross section", RuntimeWarning
            )
            return CrossSectionData(total=0, inelastic=0, elastic=0)

        # Calculate hadron-nucleon cross section using kinematic values
        if not self._cascade.sigmaSetuphN(int(kin.p1), kin.beams[0], kin.m1):
            warnings.warn(
                "Could not calculate hadron-nucleon cross section", RuntimeWarning
            )
            return CrossSectionData(total=0, inelastic=0, elastic=0)

        # Get hadron-nucleus cross section
        sigma_hA = self._cascade.sigmahA(self.p2.A)

        # For now, assume all cross section is inelastic
        # (could be refined with more detailed cross section info)
        return CrossSectionData(
            total=sigma_hA,
            inelastic=sigma_hA,
            elastic=0,  # Cascade mode focuses on inelastic processes
        )

    def _generate(self):
        """Generate a hadron-nucleus collision event."""
        # Check if cascade is initialized
        if not hasattr(self, "_cascade") or self._cascade is None:
            return False
        kin = self.kinematics
        # Set up for collision using kinematic values
        if not self._cascade.sigmaSetuphN(int(kin.p1), kin.beams[0], kin.m1):
            return False

        # Generate collision at origin (vertex can be set later)
        vertex = np.array([0, 0, 0, 0])  # (x, y, z, t)
        self._current_event = self._cascade.nextColl(kin.p2.Z, kin.p2.A, vertex)
        return self._current_event.size > 0

    def stat(self):
        """Print statistics."""
        if hasattr(self, "_cascade") and self._cascade is not None:
            self._cascade.stat()

    @property
    def max_energy(self):
        """Maximum energy for cascade initialization."""
        return self._max_energy

    @property
    def current_event(self):
        """Get the current event."""
        if self._current_event is None:
            raise RuntimeError("No current event available")
        return self._current_event
