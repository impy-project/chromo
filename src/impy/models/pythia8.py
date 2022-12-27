from impy.common import MCRun, MCEvent, CrossSectionData
from impy.util import _cached_data_dir, name2pdg
from os import environ
import numpy as np
from impy.kinematics import EventFrame
from impy.constants import standard_projectiles
from particle import literals as lp


class PYTHIA8Event(MCEvent):
    """Wrapper class around HEPEVT particle stack."""

    _hepevt = "hepevt"

    def _charge_init(self, npart):
        return self._lib.charge_from_pid(
            self._lib.pythia.particleData, self._lib.hepevt.idhep[:npart]
        )

    def _get_impact_parameter(self):
        hi = self._lib.pythia.info.hiInfo
        if hi is None:
            return np.nan
        return hi.b

    def _get_n_wounded(self):
        hi = self._lib.pythia.info.hiInfo
        if hi is None:
            return (0, 0)
        return hi.nPartProj, hi.nPartTarg


class Pythia8(MCRun):
    _name = "Pythia"
    _version = "8.308"
    _library_name = "_pythia8"
    _event_class = PYTHIA8Event
    _frame = EventFrame.CENTER_OF_MASS
    _projectiles = standard_projectiles | {lp.photon.pdgid}
    # Nuclei are supported in principle, but generation is very slow.
    # Support for nuclei can be added with ParticleData.addParticle.
    _targets = standard_projectiles | {
        name2pdg(x)
        for x in ("He4", "Li6", "C12", "O16", "Cu63", "Xe129", "Au197", "Pb208")
    }
    _restartable = True
    _data_url = (
        "https://github.com/impy-project/impy"
        + "/releases/download/zipped_data_v1.0/Pythia8_v002.zip"
    )

    def __init__(self, evt_kin, *, seed=None):
        super().__init__(seed)

        self._lib.hepevt = self._lib.Hepevt()

        datdir = _cached_data_dir(self._data_url) + "xmldoc"

        # Must delete PYTHIA8DATA from environ if it exists, since it overrides
        # our argument here. When you install Pythia8 with conda, it sets
        # PYTHIA8DATA. If that version does not match this version, readString()
        # or init() may fail.
        if "PYTHIA8DATA" in environ:
            del environ["PYTHIA8DATA"]
        self._lib.pythia = self._lib.Pythia(datdir, True)

        # must come last
        self.kinematics = evt_kin
        self._set_final_state_particles()

    def _cross_section(self):
        st = self._lib.pythia.info.sigmaTot
        return CrossSectionData(
            st.sigmaTot,
            st.sigmaTot - st.sigmaEl,
            st.sigmaEl,
            st.sigmaAX,
            st.sigmaXB,
            st.sigmaXX,
            st.sigmaAXB,
        )

    def _set_kinematics(self, kin):
        pythia = self._lib.pythia
        pythia.settings.resetAll()

        config = [
            "Random:setSeed = on",
            f"Random:seed = {self._seed}",
            # use center-of-mass frame
            "Beams:frameType = 1",
            "SoftQCD:inelastic = on",
            # reduce verbosity
            "Print:quiet = on",
            "Next:numberCount = 0",  # do not print progress
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

        for line in config:
            if not pythia.readString(line):
                raise RuntimeError(f"readString({line!r}) failed")

        # calling init several times is allowed
        if not pythia.init():
            raise RuntimeError("Pythia8 initialization failed")

    def _set_stable(self, pdgid, stable):
        self._lib.pythia.particleData.mayDecay(pdgid, not stable)

    def _generate(self):
        success = self._lib.pythia.next()
        if success:
            # We copy over the event record from Pythia's internal buffer to our
            # Hepevt-like object. This is not efficient, but easier to
            # implement. The time needed to copy the record is small compared to
            # the time needed to generate the event. If this turns out to be a
            # bottleneck, we need to make Hepevt a view of the internal record.
            self._lib.hepevt.fill(self._lib.pythia.event)
        return success
