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


class EposLHC(MCRun):
    """Implements all abstract attributes of MCRun for the
    EPOS-LHC series of event generators."""

    _name = "EPOS"
    _version = "LHC"
    _event_class = EPOSEvent
    _library_name = "_eposlhc"
    _output_frame = "center-of-mass"

    def __init__(self, event_kinematics, seed="random", logfname=None):
        from os import path

        super().__init__(seed, logfname)

        k = event_kinematics

        epos_conf = impy_config["epos"]

        info(1, "First initialization")
        self._lib.aaset(0)

        if impy_config["user_frame"] == "center-of-mass":
            iframe = 1
            self._output_frame = "center-of-mass"
        elif impy_config["user_frame"] == "laboratory":
            iframe = 2
            self._output_frame = "laboratory"

        datdir = path.join(base_path, epos_conf["datdir"])
        self._lib.initializeepos(
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
        self._set_final_state_particles()
        self._lib.charge_vect = np.vectorize(self._lib.getcharge, otypes=[np.int32])
        self._set_event_kinematics(event_kinematics)
        # Turn on particle history
        self._lib.othe1.istmax = 1

    def sigma_inel(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        return self._lib.xsection()[1]

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
        info(5, "Setting event kinematics")
        k = event_kinematics
        self._curr_event_kin = k
        self._lib.initeposevt(*self._epos_tup())

    def _attach_log(self, fname=None):
        """Routes the output to a file or the stdout."""
        fname = impy_config["output_log"] if fname is None else fname
        if fname == "stdout":
            lun = 6
            info(5, "Output is routed to stdout.")
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, "Output is routed to", fname, "via LUN", lun)

        self._lun = lun

    def _set_stable(self, pdgid, stable):
        if stable:
            self._lib.setstable(pdgid)
        else:
            self._lib.setunstable(pdgid)

    def _generate_event(self):
        self._lib.aepos(-1)
        self._lib.afinal()
        self._lib.hepmcstore()
        return False
