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

    def __init__(self, lib, event_kinematics, event_frame):
        # HEPEVT (style) common block
        evt = lib.hepevt

        # Save selector for implementation of on-demand properties
        px, py, pz, en, m = evt.phep
        vx, vy, vz, vt = evt.vhep

        self.charge_vec = None

        MCEvent.__init__(
            self,
            lib=lib,
            event_kinematics=event_kinematics,
            event_frame=event_frame,
            nevent=evt.nevhep,
            npart=evt.nhep,
            p_ids=evt.idhep,
            status=evt.isthep,
            px=px,
            py=py,
            pz=pz,
            en=en,
            m=m,
            vx=vx,
            vy=vy,
            vz=vz,
            vt=vt,
            pem_arr=evt.phep,
            vt_arr=evt.vhep,
        )

    def filter_final_state(self):
        self.selection = np.where(self.status == 1)
        self._apply_slicing()

    def filter_final_state_charged(self):
        self.selection = np.where((self.status == 1) & (self.charge != 0))
        self._apply_slicing()

    @property
    def parents(self):
        MCEvent.parents(self)
        return self.lib.hepevt.jmohep

    @property
    def children(self):
        MCEvent.children(self)
        return self.lib.hepevt.jdahep

    @property
    def _charge_init(self):
        if self.charge_vec is None:
            self.charge_vec = np.asarray([
                self.lib.pychge(self.lib.pyjets.k[i, 1]) / 3 for i in range(self.npart)
            ])
        return self.charge_vec[self.selection]


class PYTHIA6Run(MCRun):
    """Implements all abstract attributes of MCRun for the
    EPOS-LHC series of event generators."""

    def sigma_inel(self, *args, **kwargs):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        return self.lib.pyint7.sigt[0, 0, 5]

    def sigma_inel_air(self, **kwargs):
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
        self.lib.pyinit("CMS", self.p1_type, self.p2_type, self.ecm)
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

        self.lib.pydat1.mstu[10] = lun

    def init_generator(self, event_kinematics, seed="random", logfname=None):
        # from random import randint
        from impy.constants import sec2cm

        self._abort_if_already_initialized()

        # if seed == "random":
        #     seed = randint(1000000, 10000000)
        #     sseed = str(seed)
        #     self.lib.pydatr.mrpy[:4] = (
        #         int(sseed[0:2]),
        #         int(sseed[2:4]),
        #         int(sseed[4:6]),
        #         int(sseed[6:]),
        #     )
        # else:
        #     seed = int(seed)
        info(5, "Using seed:", seed)

        self.attach_log(fname=logfname)

        if impy_config["pythia6"]["new_mpi"]:
            # Latest Pythia 6 is tune 383
            self.lib.pytune(383)
            self.event_call = self.lib.pyevnw
        else:
            self.event_call = self.lib.pyevnt

        # self.mstp[51]

        # self.lib.pysubs.msel = 2
        self._set_event_kinematics(event_kinematics)

        # Set default stable
        self._define_default_fs_particles()
        # Set PYTHIA decay flags to follow all changes to MDCY
        self.lib.pydat1.mstj[21 - 1] = 1
        self.lib.pydat1.mstj[22 - 1] = 2
        # # Set ctau threshold in PYTHIA for the default stable list
        self.lib.pydat1.parj[70] = impy_config["tau_stable"] * sec2cm * 10.0  # mm

    def set_stable(self, pdgid, stable=True):
        kc = self.lib.pycomp(pdgid)
        if stable:
            self.lib.pydat3.mdcy[kc - 1, 0] = 0
            info(5, "defining", pdgid, "as stable particle")
        else:
            self.lib.pydat3.mdcy[kc - 1, 0] = 1
            info(5, pdgid, "allowed to decay")

    def generate_event(self):
        self.event_call()
        self.lib.pyhepc(1)
        return 0


class Pythia6(PYTHIA6Run):
    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.definitions import interaction_model_by_tag as models_dict

        interaction_model_def = models_dict["PYTHIA6"]
        super().__init__(interaction_model_def)
        self.init_generator(event_kinematics, seed, logfname)
