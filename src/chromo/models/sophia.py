import numpy as np
from chromo.common import MCRun, MCEvent, CrossSectionData
from chromo.constants import nucleon_mass
from chromo.constants import microbarn
from chromo.kinematics import EventFrame
from particle import literals as lp

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
    _library_name = "_sophia"
    _event_class = SophiaEvent
    _frame = EventFrame.FIXED_TARGET
    _projectiles = {lp.photon.pdgid}
    _targets = {lp.p.pdgid, lp.n.pdgid}
    _ecm_min = 0

    def _once(self, keep_decayed_particles=True):
        from chromo import debug_level

        self._lib.s_plist.ideb = debug_level
        # Keep decayed particles in the history:
        self._lib.eg_io.keepdc = keep_decayed_particles

    def _cross_section(self, kin):
        self._set_kinematics(kin)
        # code=3 for inelastic cross-section
        # TODO fill more cross-sections
        inel = (
            self._lib.crossection(self._energy_of_photon, 3, self._nucleon_code)
            * microbarn
        )
        return CrossSectionData(inelastic=inel)

    def _set_kinematics(self, kin):
        # Here we consider laboratory frame where photon moves along z axis
        # and nucleon is at rest. The angle is counted from z axis.
        # However, because of the definitions in "eventgen" subroutine of
        # SOPHIA code (line "P_gam(3) = -EPS*COS(theta*pi/180.D0)")
        # this angle should be 180 for photon moving along z
        # (and 0 for photon moving in direction opposite to z)
        self._nucleon_code = self._lib.icon_pdg_sib(kin.p2)
        self._angle_between_nucleon_and_photon = 180
        self._energy_of_photon = kin.elab
        self._energy_of_nucleon = np.float32(nucleon_mass)
        # setting parameters for cross-section
        self._lib.initial(self._nucleon_code)

    def _set_stable(self, pdgid, stable):
        sid = abs(self._lib.icon_pdg_sib(pdgid)) - 1
        idb = self._lib.s_csydec.idb
        if sid < 0 or sid > idb.size:
            raise ValueError(f"{pdgid} is unknown")
        if stable:
            idb[sid] = -abs(idb[sid])
        else:
            idb[sid] = abs(idb[sid])

    def _generate(self):
        # Generate event (the final particles and their parameters)
        # by underlying Fortran library
        self._lib.interaction_type_code = self._lib.eventgen(
            self._nucleon_code,
            self._energy_of_nucleon,
            self._energy_of_photon,
            self._angle_between_nucleon_and_photon,
        )

        # prepare hepevt common block
        self._lib.toevt()
        return True
