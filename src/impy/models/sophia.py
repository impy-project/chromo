"""
Created on 17.03.2014

@author: afedynitch
"""

import numpy as np
from impy.common import MCRun, MCEvent, impy_config
from impy.util import info

sophia_interaction_types = [
    "multipion production (fragmentation)",
    "diffractive scattering: N \u03B3 \u2192 N \u03C1",
    "direct pion production: N \u03B3 \u2192 \u0394 \u03C0",
    "direct pion production: N \u03B3 \u2192 \u0394 \u03C0",
    "diffractive scattering: N \u03B3 \u2192 N \u03C9",
    "fragmentation in resonance region",
    "excitation/decay of resonance",
]


class SophiaEvent(MCEvent):
    """Wrapper class around Sophia code"""

    def __init__(self, lib, event_kinematics, event_frame):
        evt = lib.hepevt

        # Save selector for implementation of on-demand properties
        px, py, pz, en, m = evt.phep
        # verticies (x, y, z, t) - all zeros for sophia
        vx, vy, vz, vt = evt.vhep

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
    def interaction_type(self):
        return sophia_interaction_types[self.lib.interaction_type_code]

    @property
    def _charge_init(self):
        return self.lib.schg.ichg[self.selection]

    @property
    def decayed_parent(self):
        """Returns the array of zero-based indices
        of the decayed parent particles.
        Index -1 means that there is no parent particle.
        It throw an exception (via MCEvent.parents)
        if selection is applied
        """
        MCEvent.parents(self)
        return self.lib.schg.iparnt[0 : self.npart]

    @property
    def parents(self):
        """In SOPHIA parents are difficult to obtain. This function
        returns a zeroed array of the correct shape.
        """
        MCEvent.parents(self)
        return self.lib.hepevt.jmohep[:, 0 : self.npart]

    @property
    def children(self):
        """In SOPHIA daughters are difficult to obtain. This function
        returns a zeroed array of the correct shape.
        """
        MCEvent.children(self)
        return self.lib.hepevt.jdahep[:, 0 : self.npart]


class SophiaRun(MCRun):
    """Implements all abstract attributes of MCRun for the
    Sophia event generator.
    """

    def sigma_inel(self, *args, **kwargs):
        """Inelastic cross section according to current
        event setup (energy, projectile, target).
        Currently it returns total cross section.
        """
        k = self._curr_event_kin

        if k.p1pdg != 22:
            info(0, "The first particle must be a photon, but particle pdg = ", k.p1pdg)
            raise Exception("Input error")

        if k.p2pdg not in [2212, 2112]:
            info(
                0,
                "The second particle must be a proton or neutron, but particle pdg = ",
                k.p2pdg,
            )
            raise Exception("Input error")

        # self.energy_of_photon is the energy in lab frame
        # where nucleon is at rest and photon is moving
        total_crossection_id = 3  # 3 is for total crossection
        # cross section in micro barn
        return self.lib.crossection(
            self.energy_of_photon, total_crossection_id, self.nucleon_code_number
        )

    def sigma_inel_air(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        raise Exception("SophiaRun.sigma_inel_air has no implementation")

    def _set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""

        info(5, "Setting event kinematics.")
        info(10, event_kinematics)
        k = event_kinematics

        if k.p1pdg != 22:
            info(0, "The first particle must be a photon, but particle pdg = ", k.p1pdg)
            raise Exception("Input error")

        if k.p2pdg not in [2212, 2112]:
            info(
                0,
                "The second particle must be a proton or neutron, but particle pdg = ",
                k.p2pdg,
            )
            raise Exception("Input error")

        self.nucleon_code_number = self.lib.icon_pdg_sib(k.p2pdg)
        self.energy_of_nucleon = k.pmass2
        self.energy_of_photon = k.elab
        # Here we consider laboratory frame where photon moves along z axis
        # and nucleon is at rest. The angle is counted from z axis.
        # However, because of the definitions in "eventgen" subroutine of
        # SOPHIA code (line "P_gam(3) = -EPS*COS(theta*pi/180.D0)")
        # this angle should be 180 for photon moving along z
        # (and 0 for photon moving in direction opposite to z)
        self.angle_between_nucleon_and_photon = 180
        self._curr_event_kin = event_kinematics

    def attach_log(self, fname=None):
        """Routes the output to a file or the stdout."""
        fname = impy_config["output_log"] if fname is None else fname
        if fname == "stdout":
            # self.lib.s_debug.lun = 6
            info(5, "Output is routed to stdout.")
        else:
            lun = self._attach_fortran_logfile(fname)
            # self.lib.s_debug.lun = lun
            info(5, "Output is routed to", fname, "via LUN", lun)

    def init_generator(self, event_kinematics, seed="random", logfname=None):
        from random import randint

        self._abort_if_already_initialized()

        if seed == "random":
            seed = randint(1000000, 10000000)
        else:
            seed = int(seed)
        info(5, "Using seed:", seed)

        self.lib.s_plist.ideb = impy_config["sophia"]["debug_level"]

        self._set_event_kinematics(event_kinematics)
        self.attach_log(fname=logfname)

        self.lib.init_rmmard(int(seed))  # setting random number generator seed
        self.lib.initial(
            self.nucleon_code_number
        )  # setting parameters for cross-section

        # Keep decayed particles in the history:
        self.lib.eg_io.keepdc = impy_config["sophia"]["keep_decayed_particles"]
        self._define_default_fs_particles()

    def set_stable(self, pdgid, stable=True):

        # Do not use global stable_list
        if not impy_config["sophia"]["use_stable_list"]:
            return

        sid = abs(self.lib.icon_pdg_sib(pdgid))
        if abs(pdgid) == 311:
            info(1, "Ignores K0. Use K0L/S 130/310 in final state definition.")
            return
        idb = self.lib.s_csydec.idb
        if sid == 0 or sid > idb.size - 1:
            return
        if stable:
            info(
                5, "defining as stable particle pdgid/sid = {0}/{1}".format(pdgid, sid)
            )
            idb[sid - 1] = -np.abs(idb[sid - 1])
        else:
            info(5, "pdgid/sid = {0}/{1} allowed to decay".format(pdgid, sid))
            idb[sid - 1] = np.abs(idb[sid - 1])

    def generate_event(self):
        # Generate event (the final particles and their parameters)
        # by underlying Fortran library
        self.interaction_type_code = self.lib.eventgen(
            self.nucleon_code_number,
            self.energy_of_nucleon,
            self.energy_of_photon,
            self.angle_between_nucleon_and_photon,
        )

        # pass interaction type code to MCEvent
        # via additional attribute in lib object:
        setattr(self.lib, "interaction_type_code", self.interaction_type_code)
        # prepare hepevt common block
        self.lib.toevt()
        return 0  # No rejection is implemented so far


class Sophia20(SophiaRun):
    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.definitions import interaction_model_by_tag as models_dict

        interaction_model_def = models_dict["SOPHIA20"]
        super().__init__(interaction_model_def)
        self.init_generator(event_kinematics, seed, logfname)
