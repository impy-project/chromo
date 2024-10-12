import numpy as np
from chromo.common import MCRun, MCEvent, CrossSectionData
from chromo.constants import nucleon_mass
from chromo.constants import microbarn
from chromo.kinematics import EventFrame
from particle import literals as lp
import warnings

sophia_interaction_types = [
    "multipion production (fragmentation)",
    "diffractive scattering: N \u03B3 \u2192 N \u03C1",
    "direct pion production: N \u03B3 \u2192 \u0394 \u03C0",
    "direct pion production: N \u03B3 \u2192 \u0394 \u03C0",
    "diffractive scattering: N \u03B3 \u2192 N \u03C9",
    "fragmentation in resonance region",
    "excitation/decay of resonance",
]


_sophia_unstable_pids = set(
    [
        -13,
        13,
        111,
        113,
        130,
        -211,
        211,
        -213,
        213,
        221,
        223,
        310,
        -313,
        313,
        -321,
        321,
        -323,
        323,
        331,
        333,
        -1114,
        1114,
        -2112,
        2112,
        -2114,
        2114,
        -2214,
        2214,
        -2224,
        2224,
        -3112,
        3112,
        -3114,
        3114,
        -3122,
        3122,
        -3212,
        3212,
        -3214,
        3214,
        -3222,
        3222,
        -3224,
        3224,
        -3312,
        3312,
        -3314,
        3314,
        -3322,
        3322,
        -3324,
        3324,
        -3334,
        3334,
    ]
)


class SophiaEvent(MCEvent):
    """Wrapper class around Sophia code"""

    @property
    def interaction_type(self):
        return sophia_interaction_types[self._lib.interaction_type_code]

    def _charge_init(self, npart):
        return self._lib.schg.ichg[:npart]

    def _repair_initial_beam(self):
        self._prepend_initial_beam()
        # Repair history
        # Make second mother = -1
        self.mothers[:, 1] = -1
        # Attach to initial beam particles
        self.mothers[(self.mothers == [1, -1]).all(axis=1)] = [0, 1]
        # Daughters are not defined
        self.daughters[:] = [-1, -1]


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
    _unstable_pids = _sophia_unstable_pids
    _ecm_min = 0

    def __init__(self, kinematics, *, seed=None, keep_decayed_particles=True):
        import chromo

        super().__init__(seed)

        self._lib.s_plist.ideb = chromo.debug_level
        # Keep decayed particles in the history:
        self._lib.eg_io.keepdc = keep_decayed_particles

        self.kinematics = kinematics
        self._set_final_state_particles()

    def _cross_section(self, kin=None, max_info=False):
        # code=3 for inelastic cross-section
        # TODO fill more cross-sections
        inel = (
            self._lib.crossection(self._energy_of_photon, 3, self._nucleon_code)
            * microbarn
        )
        return CrossSectionData(inelastic=inel)

    def _set_kinematics(self, evt_kin):
        # Here we consider laboratory frame where photon moves along z axis
        # and nucleon is at rest. The angle is counted from z axis.
        # However, because of the definitions in "eventgen" subroutine of
        # SOPHIA code (line "P_gam(3) = -EPS*COS(theta*pi/180.D0)")
        # this angle should be 180 for photon moving along z
        # (and 0 for photon moving in direction opposite to z)
        self._nucleon_code = self._lib.icon_pdg_sib(evt_kin.p2)
        self._angle_between_nucleon_and_photon = 180
        self._energy_of_photon = evt_kin.elab
        self._energy_of_nucleon = np.float32(nucleon_mass)
        # setting parameters for cross-section
        self._lib.initial(self._nucleon_code)

    def _set_stable(self, pdgid, stable):
        if pdgid not in self._unstable_pids:
            return

        sid = abs(self._lib.icon_pdg_sib(pdgid)) - 1
        idb = self._lib.s_csydec.idb
        if sid < 0 or sid > idb.size:
            warnings.warn(f"{pdgid} unknown to Sophia", RuntimeWarning)
            return

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
