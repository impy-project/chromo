"""
Created on 19.01.2015

@author: afedynitch
"""

from ..common import MCRun, MCEvent
from ..util import info, AZ2pdg
from impy import impy_config, base_path
from pathlib import Path
import os


class PYTHIA8Event(MCEvent):
    """Wrapper class around HEPEVT particle stack."""

    _hepevt = "hepevt"
    # We actually have daughter info, but since it is redundant, we do not copy
    # it over from the Pythia event to safe time. The time safed is probably
    # negligible though, so could be enabled if really necessary.
    _jdahep = None

    def _charge_init(self, npart):
        return self._lib.charge_from_pid(
            self._lib.pythia.particleData, self._lib.hepevt.idhep[:npart]
        )

    # Nuclear collision parameters
    @property
    def impact_parameter(self):
        """Returns impact parameter for nuclear collisions."""
        return self._lib.pythia.info.hiinfo.b

    @property
    def n_wounded_A(self):
        """Number of wounded nucleons side A"""
        return self._lib.pythia.info.hiinfo.nPartProj

    @property
    def n_wounded_B(self):
        """Number of wounded nucleons (target) side B"""
        return self._lib.pythia.info.hiinfo.nPartTarg

    @property
    def n_wounded(self):
        """Number of total wounded nucleons"""
        return (
            self._lib.pythia.info.hiinfo.nPartProj
            + self._lib.pythia.info.hiinfo.nPartTarg
        )


class Pythia8(MCRun):
    _name = "Pythia"
    _version = "8.307"
    _library_name = "_pythia8"
    _event_class = PYTHIA8Event
    _output_frame = "center-of-mass"

    def __init__(self, event_kinematics, seed=None, logfname=None):
        # store stable settings until Pythia instance is created
        self._stable = {}

        super().__init__(seed, logfname)

        self._lib.hepevt = self._lib.Hepevt()
        self._set_event_kinematics(event_kinematics)
        self._set_final_state_particles()

    def sigma_inel(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        # Cross section and energy (in mb and GeV)
        return self._pythia.info.sigmaGen()

    def _set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""

        info(5, "Setting event kinematics")

        k = event_kinematics
        self._curr_event_kin = k

        datdir = Path(base_path) / "iamdata" / "Pythia8" / "xmldoc"
        assert datdir.exists(), f"{datdir} does not exist"

        # Must delete PYTHIA8DATA from environ if it exists, since it overrides
        # our argument here. When you install Pythia8 with conda, it sets
        # PYTHIA8DATA. If that version does not match this version, readString()
        # or init() may fail.
        if "PYTHIA8DATA" in os.environ:
            del os.environ["PYTHIA8DATA"]
        pythia = self._lib.pythia = self._lib.Pythia(str(datdir), True)

        config = [
            "Random:setSeed = on",
            f"Random:seed = {self._seed}",
            # Specify energy in center of mass
            "Beams:frameType = 1",
            "SoftQCD:inelastic = on",
        ]
        # Add more options from config file
        config += impy_config["pythia8"]["options"]

        if k.p1_is_nucleus or k.p2_is_nucleus:
            config.append("HeavyIon:SigFitNGen = 0")
            config.append("HeavyIon:SigFitDefPar = 10.79,1.75,0.30,0.0,0.0,0.0,0.0,0.0")

        if k.p1_is_nucleus:
            k.p1pdg = AZ2pdg(k.A1, k.A2)

        if k.p2_is_nucleus:
            k.p2pdg = AZ2pdg(k.A2, k.Z2)

        # HD: is this necessary?
        for pid, a, z in ((k.p1pdg, k.A1, k.Z1), (k.p2pdg, k.A2, k.Z2)):
            if not pythia.particleData.isParticle(pid):
                # pdgid, p name, ap name, spin, 3*charge, color, mass
                pythia.particleData.addParticle(
                    pid,
                    f"{a * 100 + z}",
                    f"{a * 100 + z}bar",
                    1,
                    3 * int(z),
                    0,
                    float(a),
                )

        config += [
            f"Beams:idA = {k.p1pdg}",
            f"Beams:idB = {k.p2pdg}",
            f"Beams:eCM = {k.ecm}",
        ]

        for line in config:
            if not pythia.readString(line):
                raise RuntimeError(f"readString('{line}') failed")

        if not pythia.init():
            raise RuntimeError("Pythia8 initialization failed")

    def _attach_log(self, fname):
        # not implemented
        pass

    def _set_stable(self, pdgid, stable):
        self._lib.pythia.particleData.mayDecay(pdgid, not stable)

    def _generate_event(self):
        success = self._lib.pythia.next()
        if success:
            # We copy over the event record from Pythia's internal buffer to our
            # Hepevt-like object. This is not efficient, but easier to
            # implement. The time needed to copy the record is small compared to
            # the time needed to generate the event. If this turns out to be a
            # bottleneck, we need to make Hepevt a view of the interval record.
            self._lib.hepevt.fill(self._lib.pythia.event)
        return not success
