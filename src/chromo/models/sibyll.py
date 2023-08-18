from chromo.common import MCRun, MCEvent, CrossSectionData
from chromo.util import info, Nuclei
from chromo.kinematics import EventFrame
from particle import literals as lp
import warnings
import numpy as np


sibyll_decaying_pids = [
    -13,
    13,
    111,
    113,
    -211,
    211,
    -213,
    213,
    221,
    223,
    130,
    310,
    -321,
    321,
    331,
    333,
    -411,
    411,
    -413,
    413,
    -421,
    421,
    -423,
    423,
    -431,
    431,
    441,
    443,
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
    -4112,
    4112,
    -4114,
    4114,
    -4122,
    4122,
    -4132,
    4132,
    -4212,
    4212,
    -4214,
    4214,
    -4222,
    4222,
    -4224,
    4224,
    -4232,
    4232,
    -4314,
    4314,
    -4324,
    4324,
    -4332,
    4332,
]

sibyll21_decaying_pids = sibyll_decaying_pids + [-10311, 10311, -10321, 10321]
# sibyll21_decaying_pids = sibyll_decaying_pids
sibyll21_decaying_pids = set(sibyll21_decaying_pids)

sibyll23_decaying_pids = sibyll_decaying_pids + [
    -15,
    15,
    -313,
    313,
    -323,
    323,
]
sibyll23_decaying_pids = set(sibyll23_decaying_pids)


class SibyllEvent(MCEvent):
    """Wrapper class around SIBYLL 2.1 & 2.3 particle stack."""

    _jdahep = None  # no child info

    def _charge_init(self, npart):
        return self._lib.schg.ichg[:npart]

    def _get_impact_parameter(self):
        return self._lib.cnucms.b

    def _get_n_wounded(self):
        return self._lib.cnucms.na, self._lib.cnucms.nb

    @property
    def n_NN_interactions(self):
        """Number of inelastic nucleon-nucleon interactions"""
        return self._lib.cnucms.ni


class SIBYLLRun(MCRun):
    """Implements all abstract attributes of MCRun for the
    SIBYLL 2.1, 2.3 and 2.3c event generators."""

    _name = "SIBYLL"
    _event_class = SibyllEvent
    _frame = EventFrame.CENTER_OF_MASS
    _targets = Nuclei(a_max=20)
    _cross_section_projectiles = {
        p.pdgid: sib_id
        for p, sib_id in (
            (lp.p, 1),
            (lp.n, 1),
            (lp.pi_plus, 2),
            (lp.K_plus, 3),
            (lp.K_S_0, 3),
            (lp.K_L_0, 3),
        )
    }

    def __init__(self, evt_kin, *, seed=None):
        super().__init__(seed)

        # setup logging
        import chromo

        lun = 6  # stdout
        self._lib.s_debug.lun = lun
        self._lib.s_debug.ndebug = chromo.debug_level

        self._lib.sibini()
        self._lib.pdg_ini()

        # This calls _set_event_kinematics which uses self._lib.isib_pdg2pid
        # which works only after an initialization call to self._lib.pdg_ini()
        self.kinematics = evt_kin

        super()._set_final_state_particles()

    def _cross_section(self, kin=None):
        kin = self.kinematics if kin is None else kin
        if kin.p2.A > 1:
            # TODO figure out what this returns exactly:
            # self._lib.sib_sigma_hnuc
            warnings.warn(
                f"cross-section for nuclear targets not yet supported in {self.label}",
                RuntimeWarning,
            )
            return CrossSectionData()

        sib_id = self._cross_section_projectiles[abs(kin.p1)]

        tot, el, inel, diff, _, _ = self._lib.sib_sigma_hp(sib_id, kin.ecm)
        return CrossSectionData(
            total=tot,
            elastic=el,
            inelastic=inel,
            diffractive_xb=diff[0],
            diffractive_ax=diff[1],
            diffractive_xx=diff[2],
            diffractive_axb=0,
        )

    def sigma_inel_air(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        kin = self.kinematics
        sib_id = self._cross_section_projectiles[abs(kin.p1)]
        sigma = self._lib.sib_sigma_hair(sib_id, self._ecm)
        if isinstance(sigma, tuple):
            return sigma[0]
        return sigma

    def _set_kinematics(self, kin):
        self._production_id = self._lib.isib_pdg2pid(kin.p1)
        assert self._production_id != 0

    def _set_stable(self, pdgid, stable):
        if pdgid not in self._decaying_pids:
            return

        sid = abs(self._lib.isib_pdg2pid(pdgid))
        if abs(pdgid) == 311:
            info(1, "Ignores K0. Using K0L/S instead")
            self.set_stable(130, stable)
            self.set_stable(310, stable)
            return
        idb = self._lib.s_csydec.idb
        if sid == 0 or sid > idb.size - 1:
            return
        if stable:
            idb[sid - 1] = -abs(idb[sid - 1])
        else:
            idb[sid - 1] = abs(idb[sid - 1])

    def _generate(self):
        kin = self.kinematics
        self._lib.sibyll(self._production_id, kin.p2.A, kin.ecm)
        self._lib.decsib()
        self._lib.sibhep()
        return True


class Sibyll21(SIBYLLRun):
    _version = "2.1"
    _library_name = "_sib21"
    _decaying_pids = sibyll21_decaying_pids

    def _set_final_state_particles(self, pdgid):
        if not np.all(np.isin([-13, 13, -2112, 2112], pdgid)):
            raise ValueError(
                "Sibyll21 hangs when pdgs [-13, 13, -2112, 2112] are set as decaying.\n"
                "final_state_particles should contain [-13, 13, -2112, 2112]"
            )

        super()._set_final_state_particles(pdgid)


class Sibyll23(SIBYLLRun):
    _version = "2.3"
    _library_name = "_sib23"
    _decaying_pids = sibyll23_decaying_pids


class Sibyll23c(Sibyll23):
    _version = "2.3c"
    _library_name = "_sib23c01"


# undocumented patch version
class Sibyll23c00(Sibyll23):
    _version = "2.3c00"
    _library_name = "_sib23c00"


# identical to 2.3c
class Sibyll23c01(Sibyll23):
    _version = "2.3c01"
    _library_name = "_sib23c01"


# undocumented patch version
class Sibyll23c02(Sibyll23):
    _version = "2.3c02"
    _library_name = "_sib23c02"


# The c03 version was also in CORSIKA until 2020
class Sibyll23c03(Sibyll23):
    _version = "2.3c03"
    _library_name = "_sib23c03"


# The latest patch c04 was renamed to d, to generate less confusion
class Sibyll23d(Sibyll23):
    _version = "2.3d"
    _library_name = "_sib23d"
