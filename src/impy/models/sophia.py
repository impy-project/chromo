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

    @property
    def interaction_type(self):
        return sophia_interaction_types[self._lib.interaction_type_code]

    def _charge_init(self, npart):
        return self._lib.schg.ichg[:npart]

    @property
    def decayed_parent(self):
        """Returns the array of zero-based indices
        of the decayed parent particles.
        Index -1 means that there is no parent particle.
        It throw an exception (via MCEvent.parents)
        if selection is applied
        """
        return self._lib.schg.iparnt[: self.npart]


class Sophia20(MCRun):
    """Implements all abstract attributes of MCRun for the
    Sophia event generator.
    """

    _name = "Sophia"
    _version = "2.0"
    _event_class = SophiaEvent
    _library_name = "_sophia"
    _output_frame = "laboratory"

    def __init__(self, event_kinematics, seed="random", logfname=None):
        super().__init__(seed, logfname)

        self._lib.s_plist.ideb = impy_config["sophia"]["debug_level"]

        self._lib.init_rmmard(self._seed)

        # Keep decayed particles in the history:
        self._lib.eg_io.keepdc = impy_config["sophia"]["keep_decayed_particles"]
        self.event_kinematics = event_kinematics
        self._set_final_state_particles()

    @staticmethod
    def _validate_kinematics(k):
        if k.p1pdg != 22:
            raise ValueError(
                "Sophia accepts only 'gamma' as a projectile "
                + "but received pdg_id = {0}".format(k.p1pdg)
            )

        if k.p2pdg not in [2212, 2112]:
            raise ValueError(
                "Sophia accepts only 'proton' or 'neutron' as a target, "
                + "but received pdg_id = {0}".format(k.p2pdg)
            )

        # HD: is this necessary?
        if k.p2_is_nucleus:
            raise ValueError(
                "Sophia accepts only 'proton' or 'neutron' as a target, "
                + "but received nucleus = ({0}, {1})".format(k.A2, k.Z2)
            )

    def _sigma_inel(self, evt_kin):
        import warnings

        warnings.warn(
            "Returns total cross-section instead of inelastic", RuntimeWarning
        )

        self._validate_kinematics(evt_kin)

        nucleon_code_number = self._lib.icon_pdg_sib(evt_kin.p2pdg)
        energy_of_photon = self.event_kinematics.elab

        # self.energy_of_photon is the energy in lab frame
        # where nucleon is at rest and photon is moving
        total_crossection_id = 3  # 3 is for total crossection
        # cross section in micro barn
        return self._lib.crossection(
            energy_of_photon, total_crossection_id, nucleon_code_number
        )

    def _set_event_kinematics(self, k):
        info(5, "Setting event kinematics")

        self._validate_kinematics(k)

        self._nucleon_code_number = self._lib.icon_pdg_sib(k.p2pdg)

        # setting parameters for cross-section
        self._lib.initial(self._nucleon_code_number)

        # Here we consider laboratory frame where photon moves along z axis
        # and nucleon is at rest. The angle is counted from z axis.
        # However, because of the definitions in "eventgen" subroutine of
        # SOPHIA code (line "P_gam(3) = -EPS*COS(theta*pi/180.D0)")
        # this angle should be 180 for photon moving along z
        # (and 0 for photon moving in direction opposite to z)
        self._angle_between_nucleon_and_photon = 180

    def _attach_log(self, fname):
        """Routes the output to a file or the stdout."""
        fname = impy_config["output_log"] if fname is None else fname
        if fname == "stdout":
            # self._lib.s_debug.lun = 6
            info(5, "Output is routed to stdout.")
        else:
            lun = self._attach_fortran_logfile(fname)
            # self._lib.s_debug.lun = lun
            info(5, "Output is routed to", fname, "via LUN", lun)

    def _set_stable(self, pdgid, stable):

        # Do not use global stable_list
        if not impy_config["sophia"]["use_stable_list"]:
            return

        sid = abs(self._lib.icon_pdg_sib(pdgid))
        if abs(pdgid) == 311:
            info(1, "Ignores K0. Using K0L/S 130/310")
            self._set_stable(130, True)
            self._set_stable(310, True)
            return
        idb = self._lib.s_csydec.idb
        if sid == 0 or sid > idb.size - 1:
            return
        if stable:
            idb[sid - 1] = -np.abs(idb[sid - 1])
        else:
            idb[sid - 1] = np.abs(idb[sid - 1])

    def _generate_event(self):
        # Generate event (the final particles and their parameters)
        # by underlying Fortran library

        # fix roundoff error
        energy_of_nucleon = np.float32(self.event_kinematics.pmass2)
        energy_of_photon = self.event_kinematics.elab
        self._lib.interaction_type_code = self._lib.eventgen(
            self._nucleon_code_number,
            energy_of_nucleon,
            energy_of_photon,
            self._angle_between_nucleon_and_photon,
        )

        # prepare hepevt common block
        self._lib.toevt()
        return 0  # No rejection is implemented so far
