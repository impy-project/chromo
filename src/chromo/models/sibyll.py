from chromo.common import MCRun, MCEvent, CrossSectionData
from chromo.util import info, Nuclei
from chromo.kinematics import EventFrame
from chromo.constants import standard_projectiles
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

    def _history_zero_indexing(self):
        # Sibyll has only mothers
        self.mothers = self.mothers - 1

    def _repair_initial_beam(self):
        self._prepend_initial_beam()
        # Repair history
        self.mothers[(self.mothers == [1, 1]).all(axis=1)] = [0, 1]
        # Set [i, i] to [i, -1]
        condition = self.mothers[:, 0] == self.mothers[:, 1]
        self.mothers[condition, 1] = -1

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
    _projectiles = standard_projectiles | {
        3112,
        3122,
        3312,
        3322,
        3222,
        411,
        421,
        4232,
        431,
        4122,
        4132,
        4232,
        431,
        4332,
    }
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

        if hasattr(self, "_sstar_param"):
            self._lib.s_star.imod = self._sstar_param

        self._lib.sibini()
        self._lib.pdg_ini()

        # This calls _set_event_kinematics which uses self._lib.isib_pdg2pid
        # which works only after an initialization call to self._lib.pdg_ini()
        self.kinematics = evt_kin

        self._set_final_state_particles()

    def _cross_section(self, kin=None):
        kin = self.kinematics if kin is None else kin
        if kin.p1.A > 1:
            warnings.warn(
                f"Cross section for nuclear projectiles not supported in {self.label}",
                RuntimeWarning,
            )
            return CrossSectionData()

        sib_id = self._cross_section_projectiles[abs(kin.p1)]

        if kin.p2.A > 19:
            raise ValueError(
                f"{self.label} does not support nuclear targets heavier than 19"
            )
        if kin.p2.A > 1 and self.version == "2.1":
            raise ValueError(
                f"{self.label} 2.1 does not (yet) support nuclear targets."
            )

        if kin.p2.A > 1:
            alam = 1.0  # Not used
            icsmod = 1  # use Sibyll p-p cross section as input
            iparm = 2  # use Goulianos param. for inel. coupling param.
            self._lib.sig_had_nuc(sib_id, kin.p2.A, kin.ecm, alam, icsmod, iparm)
            nsig = self._lib.nucsig
            return CrossSectionData(
                total=float(nsig.sigt),
                prod=float(nsig.sigt - nsig.sigqe),
                quasielastic=float(nsig.sigqe),
                inelastic=float(nsig.siginel),
                diffractive_sum=float(nsig.sigqsd),
                diffractive_xb=float(nsig.sigsd),
                elastic=float(nsig.sigel),
            )

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
        sigma = self._lib.sib_sigma_hair(sib_id, kin.ecm)
        if isinstance(sigma, tuple):
            return sigma[0]
        return sigma

    def _set_kinematics(self, kin):
        self._production_id = self._lib.isib_pdg2pid(kin.p1)
        assert self._production_id != 0

    def _set_stable(self, pdgid, stable):
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
    _projectiles = standard_projectiles
    _library_name = "_sib21"


class Sibyll23(SIBYLLRun):
    _version = "2.3"
    _projectiles = standard_projectiles
    _library_name = "_sib23"


class Sibyll23c(SIBYLLRun):
    _version = "2.3c"
    _library_name = "_sib23c01"


# undocumented patch version
class Sibyll23c00(SIBYLLRun):
    _version = "2.3c00"
    _library_name = "_sib23c00"


# identical to 2.3c
class Sibyll23c01(SIBYLLRun):
    _version = "2.3c01"
    _library_name = "_sib23c01"


# undocumented patch version
class Sibyll23c02(SIBYLLRun):
    _version = "2.3c02"
    _library_name = "_sib23c02"


# The c03 version was also in CORSIKA until 2020
class Sibyll23c03(SIBYLLRun):
    _version = "2.3c03"
    _library_name = "_sib23c03"


# The latest patch c04 was renamed to d, to generate less confusion
class Sibyll23d(SIBYLLRun):
    _version = "2.3d"
    _library_name = "_sib23d"


class Sibyll23StarNoEnh(Sibyll23d):
    _version = "2.3Star-noenh"
    _library_name = "_sib23d_star"
    _sstar_param = 0


class Sibyll23StarRho(Sibyll23StarNoEnh):
    _version = "2.3Star-rho"
    _library_name = "_sib23d_star"
    _sstar_param = 1


class Sibyll23StarBar(Sibyll23StarNoEnh):
    _version = "2.3Star-bar"
    _library_name = "_sib23d_star"
    _sstar_param = 2


class Sibyll23StarStrange(Sibyll23StarNoEnh):
    _version = "2.3Star-strange"
    _library_name = "_sib23d_star"
    _sstar_param = 3


class Sibyll23StarMixed(Sibyll23StarNoEnh):
    _version = "2.3Star-mix"
    _library_name = "_sib23d_star"
    _sstar_param = 4
