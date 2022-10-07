"""
Created on 19.01.2015

@author: afedynitch
"""

import numpy as np
from impy.common import MCRun, MCEvent
from impy import impy_config
from impy.util import info, AZ2pdg


class PYTHIA8Event(MCEvent):
    """Wrapper class around HEPEVT particle stack."""

    _hepevt = "hepevt"
    _jdahep = None

    def _charge_init(self, npart):
        # TODO
        return np.zeros(npart, dtype=np.int_)

    # Nuclear collision parameters
    @property
    def impact_parameter(self):
        """Returns impact parameter for nuclear collisions."""
        return self._lib.info.hiinfo.b()

    @property
    def n_wounded_A(self):
        """Number of wounded nucleons side A"""
        return self._lib.info.hiinfo.nPartProj()

    @property
    def n_wounded_B(self):
        """Number of wounded nucleons (target) side B"""
        return self._lib.info.hiinfo.nPartTarg()

    @property
    def n_wounded(self):
        """Number of total wounded nucleons"""
        return self._lib.info.hiinfo.nPartProj() + self._lib.info.hiinfo.nPartTarg()


class PYTHIA8Run(MCRun):
    def sigma_inel(self, *args, **kwargs):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        # Cross section and energy (in mb and GeV)
        return self.lib.info.sigmaGen()

    def _set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""
        k = event_kinematics
        self._curr_event_kin = k

        # create new object
        self.lib = self.cpp_lib.Pythia("", True)
        # Replay initialization strings
        for param_string in self.save_init_strings:
            self.lib.readString(param_string)
        # Replay stable history
        for pdgid in self.stable_history:
            self.set_stable(pdgid, self.stable_history[pdgid])

        if k.p1_is_nucleus or k.p2_is_nucleus:
            self.lib.readString("HeavyIon:SigFitNGen = 0")
            self.lib.readString(
                "HeavyIon:SigFitDefPar = 10.79,1.75,0.30,0.0,0.0,0.0,0.0,0.0"
            )
        if k.p1_is_nucleus:
            k.p1pdg = AZ2pdg(k.A1, k.A2)
            # pdgid, p name, ap name, spin, 3*charge, color, mass
            self.lib.particleData.addParticle(
                k.p1pdg,
                str(k.A1 * 100 + k.Z1),
                str(k.A1 * 100 + k.Z1) + "bar",
                1,
                3 * k.Z1,
                0,
                float(k.A1),
            )
        if k.p2_is_nucleus:
            k.p2pdg = AZ2pdg(k.A2, k.Z2)
            # pdgid, p name, ap name, spin, 3*charge, color, mass
            self.lib.particleData.addParticle(
                k.p2pdg,
                str(k.A2 * 100 + k.Z2),
                str(k.A2 * 100 + k.Z2) + "bar",
                1,
                3 * k.Z2,
                0,
                float(k.A2),
            )

        self.lib.readString(f"Beams:idA = {k.p1pdg}")
        self.lib.readString(f"Beams:idB = {k.p2pdg}")
        self.lib.readString(f"Beams:eCM = {k.ecm}")
        # Set default stable
        self._define_default_fs_particles()

        self.lib.init()

        info(5, "Setting event kinematics")

    def attach_log(self, fname):
        fname = impy_config["output_log"] if fname is None else fname

    def init_generator(self, event_kinematics, seed="random", logfname=None):
        from random import randint

        self._abort_if_already_initialized()

        if seed == "random":
            seed = randint(1000000, 10000000)
        else:
            seed = int(seed)
        info(5, "Using seed:", seed)

        # Since a Pythia 8 instance is an object unlike in the case
        # of the Fortran stuff where the import of self.lib generates
        # the object, we will backup the library
        self.cpp_lib = self.lib
        self.lib = None
        # The object/instance is created each time event_kinematics is
        # set, since Pythia8 ATM does not support changing beams or
        # energies at runtime.
        # Super cool workaround but not very performant!

        self.save_init_strings = [
            "Random:setSeed = on",
            "Random:seed = " + str(seed),
            # Specify energy in center of mass
            "Beams:frameType = 1",
            # Minimum bias events
            "SoftQCD:all = on",
        ]

        # Add more options from config file
        for param_string in impy_config["pythia8"]["options"]:
            info(5, "Using Pythia 8 parameter:", param_string)
            self.save_init_strings.append(param_string)

        self.attach_log(fname=logfname)

        # Create history for particle decay settings
        # since changing energy changes resets the stable settings
        self.stable_history = {}

        self._set_event_kinematics(event_kinematics)

    def set_stable(self, pdgid, stable=True):

        may_decay = not stable
        if self.lib is not None:
            self.lib.particleData.mayDecay(pdgid, may_decay)

        info(5, pdgid, "allowed to decay: ", may_decay)

        self.stable_history[pdgid] = stable

    def generate_event(self):
        return not self.lib.next()


class Pythia8(PYTHIA8Run):
    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.definitions import interaction_model_by_tag as models_dict

        interaction_model_def = models_dict["PYTHIA8"]
        super().__init__(interaction_model_def)
        self.init_generator(event_kinematics, seed, logfname)
