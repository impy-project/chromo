"""
Created on 19.01.2015

@author: afedynitch
"""

import numpy as np
from impy.common import MCRun, MCEvent
from impy import impy_config, pdata
from impy.util import info


class PYTHIA6Event(MCEvent):
    """Wrapper class around HEPEVT particle stack."""

    def _charge_init(self, npart):
        k = self._lib.pyjets.k[:npart, 1]
        # TODO accelerate by implementing this loop in Fortran
        return np.fromiter((self._lib.pychge(ki) / 3 for ki in k), np.double)


class Pythia6(MCRun):
    """Implements all abstract attributes of MCRun for the
    EPOS-LHC series of event generators."""

    _name = "Pythia"
    _version = "6.428"
    _library_name = "_pythia6"
    _event_class = PYTHIA6Event
    _output_frame = "center-of-mass"

    def __init__(self, event_kinematics, seed=None, logfname=None):
        from impy.constants import sec2cm

        super().__init__(seed, logfname)

        self._lib.init_rmmard(self._seed)

        if impy_config["pythia6"]["new_mpi"]:
            # Latest Pythia 6 is tune 383
            self._lib.pytune(383)
            self.event_call = self._lib.pyevnw
        else:
            self.event_call = self._lib.pyevnt

        # self.mstp[51]

        # Setup pythia processes (set to custom mode)
        self._lib.pysubs.msel = 0

        # Enable minimum bias processes incl diffraction, low-pt
        # but no elastic (see p227 of hep-ph/0603175)
        for isub in [11, 12, 13, 28, 53, 68, 92, 93, 94, 95, 96]:
            self._lib.pysubs.msub[isub - 1] = 1

        self._set_event_kinematics(event_kinematics)

        # Set default stable
        self._set_final_state_particles()
        # Set PYTHIA decay flags to follow all changes to MDCY
        self._lib.pydat1.mstj[21 - 1] = 1
        self._lib.pydat1.mstj[22 - 1] = 2
        # # Set ctau threshold in PYTHIA for the default stable list
        self._lib.pydat1.parj[70] = impy_config["tau_stable"] * sec2cm * 10.0  # mm

    def sigma_inel(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        return self._lib.pyint7.sigt[0, 0, 5]

    def sigma_inel_air(self):
        """PYTHIA6 does not support nuclear targets."""
        raise Exception("PYTHIA6 does not support nuclear targets.")

    def _set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""
        k = event_kinematics
        self._curr_event_kin = k
        self.p1_type = pdata.name(k.p1pdg)
        self.p2_type = pdata.name(k.p2pdg)
        self.ecm = k.ecm
        self._lib.pyinit("CMS", self.p1_type, self.p2_type, self.ecm)
        info(5, "Setting event kinematics")

    def _attach_log(self, fname=None):
        """Routes the output to a file or the stdout."""
        fname = impy_config["output_log"] if fname is None else fname
        if fname == "stdout":
            lun = 6
            info(5, "Output is routed to stdout.")
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, "Output is routed to", fname, "via LUN", lun)

        self._lib.pydat1.mstu[10] = lun

    def _set_stable(self, pdgid, stable):
        kc = self._lib.pycomp(pdgid)
        self._lib.pydat3.mdcy[kc - 1, 0] = 0 if stable else 1

    def _generate_event(self):
        self.event_call()
        self._lib.pyhepc(1)
        return False
