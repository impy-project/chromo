from chromo.common import MCRun, EventData, CrossSectionData
from chromo.util import _cached_data_dir, name2pdg
from os import environ
import numpy as np
from chromo.kinematics import EventFrame
from chromo.constants import standard_projectiles
from particle import literals as lp
import warnings
from typing import Collection, List


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
            f"{model}-{version}: only part of the history " "available in HepMC3 event",
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
    _version = "8.308"
    _library_name = "_pythia8"
    _event_class = PYTHIA8Event
    _frame = EventFrame.CENTER_OF_MASS
    _projectiles = standard_projectiles | {lp.photon.pdgid}
    # Nuclei are supported in principle, but generation is very slow.
    # Support for more nuclei can be added with ParticleData.addParticle.
    _targets = _projectiles | {
        name2pdg(x)
        for x in ("He4", "Li6", "C12", "O16", "Cu63", "Xe129", "Au197", "Pb208")
    }
    _restartable = True
    _data_url = (
        "https://github.com/impy-project/chromo"
        + "/releases/download/zipped_data_v1.0/Pythia8_v002.zip"
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
                raise RuntimeError(f"readString({line!r}) failed")

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
    def _parse_config(config) -> List[str]:
        # convert config to lines and filter out lines that
        # we override with our settings in _set_kinematics
        result: List[str] = []
        if isinstance(config, str):
            for line in config.split("\n"):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                result.append(line)
        elif isinstance(config, Collection):
            for item in config:
                result.append(item.strip())

        ignored = ("Random:", "Beams:", "Next:")

        result2: List[str] = []
        for line in result:
            line = line.strip()
            for ig in ignored:
                if line.startswith(ig):
                    warnings.warn(f"configuration ignored: {line!r}", RuntimeWarning)
                    break
            result2.append(line)
        return result2
