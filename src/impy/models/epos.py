"""
Created on 03.05.2016

@author: afedynitch
"""

import numpy as np
from impy.common import MCRun, MCEvent
from impy import impy_config, base_path
from impy.util import info


class EPOSEvent(MCEvent):
    """Wrapper class around EPOS particle stack."""

    def __init__(self, generator):
        super().__init__(generator)
        # EPOS sets parents of beam particles to (-1, -1).
        # We change it to (0, 0)
        nbeam = np.sum(self.status == 4)
        self.parents[:nbeam] = 0

    def _charge_init(self, npart):
        return self._lib.charge_vect(self._lib.hepevt.idhep[:npart])

    # Nuclear collision parameters
    @property
    def impact_parameter(self):
        """Returns impact parameter for nuclear collisions."""
        # return self._lib.nuc3.bimp
        return self._lib.cevt.bimevt

    @property
    def n_wounded_A(self):
        """Number of wounded nucleons side A"""
        return self._lib.cevt.npjevt

    @property
    def n_wounded_B(self):
        """Number of wounded nucleons (target) side B"""
        return self._lib.cevt.ntgevt

    @property
    def n_wounded(self):
        """Number of total wounded nucleons"""
        return self._lib.cevt.npjevt + self._lib.cevt.ntgevt

    @property
    def n_spectator_A(self):
        """Number of spectator nucleons side A"""
        return self._lib.cevt.npnevt + self._lib.cevt.nppevt

    @property
    def n_spectator_B(self):
        """Number of spectator nucleons (target) side B"""
        return self._lib.cevt.ntnevt + self._lib.cevt.ntpevt


class EPOSRun(MCRun):
    """Implements all abstract attributes of MCRun for the
    EPOS-LHC series of event generators."""

    def sigma_inel(self, *args, **kwargs):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        return self.lib.xsection()[1]

    def _epos_tup(self):
        """Constructs an tuple of arguments for calls to event generator
        from given event kinematics object."""
        k = self._curr_event_kin
        info(
            20,
            "Request EPOS ARGs tuple:\n",
            (k.ecm, -1.0, k.p1pdg, k.p2pdg, k.A1, k.Z1, k.A2, k.Z2),
        )
        return (k.ecm, -1.0, k.p1pdg, k.p2pdg, k.A1, k.Z1, k.A2, k.Z2)

    def _set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""
        k = event_kinematics
        self._curr_event_kin = k
        self.lib.initeposevt(*self._epos_tup())
        info(5, "Setting event kinematics")

    def attach_log(self, fname=None):
        """Routes the output to a file or the stdout."""
        fname = impy_config["output_log"] if fname is None else fname
        if fname == "stdout":
            lun = 6
            info(5, "Output is routed to stdout.")
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, "Output is routed to", fname, "via LUN", lun)

        self._lun = lun

    def init_generator(self, event_kinematics, seed="random", logfname=None):
        from random import randint
        from os import path

        self._abort_if_already_initialized()
        k = event_kinematics
        if seed == "random":
            seed = randint(1000000, 10000000)
        else:
            seed = int(seed)
        info(5, "Using seed:", seed)

        epos_conf = impy_config["epos"]
        datdir = path.join(base_path, epos_conf["datdir"])
        self.attach_log(fname=logfname)
        info(1, "First initialization")
        self.lib.aaset(0)

        if impy_config["user_frame"] == "center-of-mass":
            iframe = 1
            self._output_frame = "center-of-mass"
        elif impy_config["user_frame"] == "laboratory":
            iframe = 2
            self._output_frame = "laboratory"

        self.lib.initializeepos(
            float(seed),
            k.ecm,
            datdir,
            len(datdir),
            iframe,
            k.p1pdg,
            k.p2pdg,
            k.A1,
            k.Z1,
            k.A2,
            k.Z2,
            impy_config["epos"]["debug_level"],
            self._lun,
        )

        # Set default stable
        self._define_default_fs_particles()
        self.lib.charge_vect = np.vectorize(self.lib.getcharge, otypes=[np.int32])
        self._set_event_kinematics(event_kinematics)
        # Turn on particle history
        self.lib.othe1.istmax = 1

    def set_stable(self, pdgid, stable=True):
        if stable:
            self.lib.setstable(pdgid)
            info(5, "defining", pdgid, "as stable particle")
        else:
            self.lib.setunstable(pdgid)
            info(5, pdgid, "allowed to decay")

    def generate_event(self):
        self.lib.aepos(-1)
        self.lib.afinal()
        self.lib.hepmcstore()
        return False


class EposLHC(EPOSRun):
    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.definitions import interaction_model_by_tag as models_dict

        interaction_model_def = models_dict["EPOSLHC"]
        super().__init__(interaction_model_def)
        self.init_generator(event_kinematics, seed, logfname)
