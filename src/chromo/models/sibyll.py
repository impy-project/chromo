import warnings

import numpy as np
from particle import literals as lp

from chromo.common import CrossSectionData, MCEvent, MCRun
from chromo.constants import standard_projectiles
from chromo.kinematics import EventFrame
from chromo.util import Nuclei, info, is_real_nucleus

_sibyll_unstable_pids = [
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

_sibyll21_unstable_pids = {*_sibyll_unstable_pids, -10311, 10311, -10321, 10321}

_sibyll23_unstable_pids = {*_sibyll_unstable_pids, -15, 15, -313, 313, -323, 323}


class SibyllEvent(MCEvent):
    """Wrapper class around SIBYLL 2.1 & 2.3 particle stack."""

    _jdahep = None  # no child info

    def _get_charge(self, npart):
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
    _targets = Nuclei(a_max=20)
    _glauber_trials = 1000  # default number of trials for Glauber model integration
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
    _projectiles = (
        standard_projectiles
        | {
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
        | Nuclei(a_max=56)
    )
    _unstable_pids = _sibyll23_unstable_pids

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

        super()._set_final_state_particles()

    def _cross_section(self, kin=None, max_info=False):
        kin = self.kinematics if kin is None else kin
        if abs(kin.p1) not in self._cross_section_projectiles and not is_real_nucleus(
            kin.p1
        ):
            warnings.warn(
                f"Cross section for {kin.p1} projectiles not supported in {self.label}",
                RuntimeWarning,
            )
            return CrossSectionData()

        if kin.p2.A > 20:
            msg = f"{self.label} does not support nuclear targets heavier than 20"
            raise ValueError(msg)
        sib_id = self._cross_section_projectiles.get(abs(kin.p1), None)

        if kin.p2.A > 1 and not is_real_nucleus(kin.p1):
            if self.version == "2.1":
                totpp, _, _, _, slopepp, rhopp = self._lib.sib_sigma_hp(sib_id, kin.ecm)
                sigt, sigel, sigqe = self._lib.glauber(kin.p2.A, totpp, slopepp, rhopp)
                return CrossSectionData(
                    total=float(sigt),
                    inelastic=float(sigt - sigel),
                    prod=float(sigt - sigqe),
                    quasielastic=float(sigqe),
                    elastic=float(sigel),
                )
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

        # nucleus-nucleon collisions
        if is_real_nucleus(kin.p1):
            # Nucleus-nucleus collisions
            self._lib.sigma_nuc_nuc(kin.p1.A, kin.p2.A, kin.ecm, self._glauber_trials)
            nsig = self._lib.nucnucsig
            return CrossSectionData(
                prod=float(nsig.sigprod),
                quasielastic=float(nsig.sigqe),
            )

        # hadron-nucleon collisions
        if sib_id is None:
            msg = f"Variable sib_id expected to be set for {kin.p1}."
            raise ValueError(msg)
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

    @property
    def glauber_trials(self):
        """Number of trials for Glauber model integration

        Default is 1000 (set at model initialisation).
        Larger number of `ntrials` reduces the fluctuations in the cross section,
        thus, making it more smooth. Smaller number of `ntrials` makes calculations of
        cross section faster.
        """
        return self._glauber_trials

    @glauber_trials.setter
    def glauber_trials(self, ntrials):
        self._glauber_trials = ntrials

    def _set_kinematics(self, kin):
        if not kin.p1.is_nucleus:
            self._production_id = self._lib.isib_pdg2pid(kin.p1)
            if self._production_id == 0:
                msg = f"Invalid _production_id: {self._production_id}. Check the input kinematics: {kin.p1}"
                raise ValueError(msg)

    def _set_stable(self, pdgid, stable):
        if pdgid not in self._unstable_pids:
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
        if kin.p1.is_nucleus:
            # Nucleus-nucleus collisions
            self._lib.sibnuc(kin.p1.A, kin.p2.A, kin.ecm)
        else:
            # hadron-nucleon collisions
            self._lib.sibyll(self._production_id, kin.p2.A, kin.ecm)
        self._lib.decsib()
        self._lib.sibhep()
        return True


class Sibyll21(SIBYLLRun):
    _version = "2.1"
    _projectiles = standard_projectiles | Nuclei(a_max=56)
    _library_name = "_sib21"
    _unstable_pids = _sibyll21_unstable_pids

    def _set_final_state_particles(self, pdgid):
        if not np.all(np.isin([-13, 13, -2112, 2112], pdgid)):
            raise ValueError(
                "Sibyll21 hangs when pdgs [-13, 13, -2112, 2112] are set as unstable.\n"
                "final_state_particles should contain [-13, 13, -2112, 2112]"
            )

        super()._set_final_state_particles(pdgid)


class Sibyll23(SIBYLLRun):
    _version = "2.3"
    _projectiles = standard_projectiles | Nuclei(a_max=56)
    _library_name = "_sib23"


class Sibyll23c(SIBYLLRun):
    _version = "2.3c"
    _library_name = "_sib23c01"


# undocumented patch version
class Sibyll23c00(Sibyll23c):
    _version = "2.3c00"
    _library_name = "_sib23c00"


# identical to 2.3c
class Sibyll23c01(Sibyll23c):
    _version = "2.3c01"
    _library_name = "_sib23c01"


# undocumented patch version
class Sibyll23c02(Sibyll23c):
    _version = "2.3c02"
    _library_name = "_sib23c02"


# The c03 version was also in CORSIKA until 2020
class Sibyll23c03(Sibyll23c):
    _version = "2.3c03"
    _library_name = "_sib23c03"


# The latest patch c04 was renamed to d, to generate less confusion
class Sibyll23d(Sibyll23c):
    _version = "2.3d"
    _library_name = "_sib23d"


class Sibyll23e(Sibyll23c):
    _version = "2.3e"
    _library_name = "_sib23e"


class Sibyll23dStarNoEnh(Sibyll23d):
    _version = "2.3dStar-noenh"
    _library_name = "_sib23d_star"
    _sstar_param = 0


class Sibyll23dStarRho(Sibyll23dStarNoEnh):
    _version = "2.3dStar-rho"
    _sstar_param = 1


class Sibyll23dStarBar(Sibyll23dStarNoEnh):
    _version = "2.3dStar-bar"
    _sstar_param = 2


class Sibyll23dStarStrange(Sibyll23dStarNoEnh):
    _version = "2.3dStar-strange"
    _sstar_param = 3


class Sibyll23dStarMixed(Sibyll23dStarNoEnh):
    _version = "2.3dStar-mix"
    _sstar_param = 4


class Sibyll23eStarNoEnh(Sibyll23e):
    _version = "2.3eStar-noenh"
    _library_name = "_sib23e_star"
    _sstar_param = 0


class Sibyll23eStarRho(Sibyll23eStarNoEnh):
    _version = "2.3eStar-rho"
    _sstar_param = 1


class Sibyll23eStarBar(Sibyll23eStarNoEnh):
    _version = "2.3eStar-bar"
    _sstar_param = 2


class Sibyll23eStarStrange(Sibyll23eStarNoEnh):
    _version = "2.3eStar-strange"
    _sstar_param = 3


class Sibyll23eStarMixed(Sibyll23eStarNoEnh):
    _version = "2.3eStar-mix"
    _sstar_param = 4
