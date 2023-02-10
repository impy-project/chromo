from chromo.common import Model as MCRun, MCEvent, CrossSectionData
from chromo.util import Nuclei
from chromo.remote_control import MCRunRemote
from chromo.kinematics import EventFrame
from particle import literals as lp
import warnings


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


class SIBYLLRun:
    """Implements all abstract attributes of MCRun for the
    SIBYLL 2.1, 2.3 and 2.3c event generators."""

    _name = "SIBYLL"
    _event_class = SibyllEvent
    _frame = EventFrame.CENTER_OF_MASS
    _targets = Nuclei(a_max=20)

    def _once(self):
        from chromo import debug_level

        # setup logging

        lun = 6  # stdout
        self._lib.s_debug.lun = lun
        self._lib.s_debug.ndebug = debug_level

        self._lib.sibini()
        self._lib.pdg_ini()

    def _cross_section(self, kin):
        self._set_kinematics(kin)
        if self._target_a > 1:
            # TODO figure out what this returns exactly:
            # self._lib.sib_sigma_hnuc
            warnings.warn(
                f"cross-section for nuclear targets not yet supported in {self.label}",
                RuntimeWarning,
            )
            return CrossSectionData()

        tot, el, inel, diff, _, _ = self._lib.sib_sigma_hp(
            self._projectile_class, self._ecm
        )
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
        sigma = self._lib.sib_sigma_hair(self._projectile_class, self._ecm)
        if isinstance(sigma, tuple):
            return sigma[0]
        return sigma

    def _set_kinematics(self, kin):
        self._projectile_id = self._lib.isib_pdg2pid(kin.p1)
        self._projectile_class = {
            lp.p.pdgid: 1,
            lp.n.pdgid: 1,
            lp.pi_plus.pdgid: 2,
            lp.K_plus.pdgid: 3,
            lp.K_S_0.pdgid: 3,
            lp.K_L_0.pdgid: 3,
        }[abs(kin.p1)]
        self._target_a = kin.p2.A
        assert self._projectile_id != 0
        assert self._target_a >= 1
        self._ecm = kin.ecm

    def _set_stable(self, pdgid, stable):
        sid = abs(self._lib.isib_pdg2pid(pdgid))
        idb = self._lib.s_csydec.idb
        if sid == 0 or sid > idb.size - 1:
            return
        if stable:
            idb[sid - 1] = -abs(idb[sid - 1])
        else:
            idb[sid - 1] = abs(idb[sid - 1])

    def _generate(self):
        self._lib.sibyll(self._projectile_id, self._target_a, self._ecm)
        self._lib.decsib()
        self._lib.sibhep()
        return True


class Sibyll21(SIBYLLRun, MCRun):
    _version = "2.1"
    _library_name = "_sib21"


# For some reason, Sibyll23 requires MCRunRemote, but the others don't
class Sibyll23(SIBYLLRun, MCRunRemote):
    _version = "2.3"
    _library_name = "_sib23"


class Sibyll23c(SIBYLLRun, MCRun):
    _version = "2.3c"
    _library_name = "_sib23c01"


# undocumented patch version
class Sibyll23c00(SIBYLLRun, MCRun):
    _version = "2.3c00"
    _library_name = "_sib23c00"


# identical to 2.3c
class Sibyll23c01(SIBYLLRun, MCRun):
    _version = "2.3c01"
    _library_name = "_sib23c01"


# undocumented patch version
class Sibyll23c02(SIBYLLRun, MCRun):
    _version = "2.3c02"
    _library_name = "_sib23c02"


# The c03 version was also in CORSIKA until 2020
class Sibyll23c03(SIBYLLRun, MCRun):
    _version = "2.3c03"
    _library_name = "_sib23c03"


# The latest patch c04 was renamed to d, to generate less confusion
class Sibyll23d(SIBYLLRun, MCRun):
    _version = "2.3d"
    _library_name = "_sib23d"
