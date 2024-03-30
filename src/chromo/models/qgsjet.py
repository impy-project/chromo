import numpy as np
from chromo.common import MCRun, MCEvent, CrossSectionData
from chromo.kinematics import EventFrame
from chromo.constants import standard_projectiles
from chromo.util import _cached_data_dir, Nuclei
from particle import literals as lp


class QGSJET1Event(MCEvent):
    """Wrapper class around QGSJet HEPEVT converter."""

    def _charge_init(self, npart):
        return self._lib.qgchg.ichg[:npart]

    @property
    def diffr_type(self):
        """Type of diffration"""
        return self._lib.jdiff.jdiff

    def _repair_initial_beam(self):
        self._prepend_initial_beam()
        # Repair history
        self.mothers[(self.mothers == [1, 1]).all(axis=1)] = [0, 1]
        # Set [i, i] to [i, -1]
        condition = self.mothers[:, 0] == self.mothers[:, 1]
        self.mothers[condition, 1] = -1
        # No daughters
        self.daughters[:] = [-1, -1]


class QGSJET2Event(QGSJET1Event):
    def _get_impact_parameter(self):
        return self._lib.qgarr7.b

    def _get_n_wounded(self):
        return self._lib.qgarr55.nwp, self._lib.qgarr55.nwt


class QGSJetRun(MCRun):
    _name = "QGSJet"
    _frame = EventFrame.FIXED_TARGET
    _projectiles = standard_projectiles | Nuclei()
    _data_url = (
        "https://github.com/impy-project/chromo"
        + "/releases/download/zipped_data_v1.0/qgsjet_v001.zip"
    )

    def __init__(self, evt_kin, *, seed=None):
        import chromo

        super().__init__(seed)

        # logging
        lun = 6  # stdout
        datdir = _cached_data_dir(self._data_url)
        self._lib.cqgsini(datdir, lun, chromo.debug_level)

        self.kinematics = evt_kin
        self._set_final_state_particles()
        self._activate_decay_handler(on=True)

    def _set_stable(self, pdgid, stable):
        # use Pythia8 instance to decay particles which QGSJet does not decay
        pass


class QGSJet1Run(QGSJetRun):
    """Implements all abstract attributes of MCRun for the
    QGSJET-01c legacy event generators."""

    _event_class = QGSJET1Event

    def _tabulated_cross_section(self, kin=None):
        # Interpolation routine for QGSJET01D cross sections from CORSIKA.
        from scipy.interpolate import UnivariateSpline

        kin = self.kinematics if kin is None else kin
        if kin.p1.A is not None and kin.p1.A > 1:
            return CrossSectionData(prod=self._lib.sectnu(kin.elab, kin.p1.A, kin.p2.A))

        A_target = kin.p2.A
        # Projectile ID-1 to access fortran indices directly
        icz = self._projectile_id - 1
        qgsgrid = 10 ** np.arange(1, 11)
        cross_section = np.zeros(10)
        wa = np.zeros(3)
        for je in range(10):
            sectn = 0.0

            ya = A_target

            ya = np.log(ya) / 1.38629 + 1.0
            ja = min(int(ya), 2)
            wa[1] = ya - ja
            wa[2] = wa[1] * (wa[1] - 1) * 0.5
            wa[0] = 1.0 - wa[1] + wa[2]
            wa[1] = wa[1] - 2.0 * wa[2]
            sectn = sum(
                [self._lib.xsect.gsect[je, icz, ja + m - 1] * wa[m] for m in range(3)]
            )
            cross_section[je] = np.exp(sectn)

        spl = UnivariateSpline(
            np.log(qgsgrid), np.log(cross_section), ext="extrapolate", s=0, k=1
        )

        prod = np.exp(spl(np.log(kin.elab)))
        return CrossSectionData(prod=prod)

    def _cross_section(self, kin=None, max_info=False):
        """
        Calculate cross section for QGSJET-01c event generator.

        Parameters
        ----------
        kin : EventKinematics, optional
        max_info : bool, optional
            Flag indicating whether to calculate the cross section for
            AA collisions using MC. Returns the cross sections for
            specific processes. Default is False and only inelastic or
            production cross section is calculated.

        Returns
        -------
        CrossSectionData
            An object containing the calculated cross section data.

        Notes
        -----
        This method calculates the cross section for the QGSJET-01c
        event generator. The cross section is calculated based on the
        provided kinematics and collision type. If `max_info is True`
        and the projectile particle is a nucleus (A > 1), the cross
        section is calculated using the Glauber model. Otherwise, the
        cross section is calculated using the tabulated cross section.
        For h-A collisions, the cross section is calculated using the
        tabulated cross section due to an issue in calculating it from
        QGSJET directly.
        """
        kin = self.kinematics if kin is None else kin

        if max_info and (kin.p1.is_nucleus and kin.p1.A > 1):
            gtot, gprod, _, gdd, gqel, gcoh = self._lib.crossc(self.glauber_trials)
            return CrossSectionData(
                total=gtot,
                prod=gprod,
                diffractive_xb=gdd,
                quasielastic=gqel,
                coherent=gcoh,
            )
        elif (kin.p1.is_nucleus and kin.p1.A > 1) or (
            kin.p2.is_nucleus and kin.p2.A > 1
        ):
            # Note: this is a workaround until the calculated h-A
            # cross sections work correctly
            return self._tabulated_cross_section(kin)

        gtot, gin, gel, gdp, gdt, gdd = self._lib.cqgshh_ha_cs()
        return CrossSectionData(
            total=gtot,
            inelastic=gin,
            elastic=gel,
            diffractive_xb=gdp,
            diffractive_ax=gdt,
            diffractive_xx=gdd,
        )

    @property
    def glauber_trials(self):
        """Number of trials for Glauber model integration

        Default is 1000 (set at model initialisation).
        Larger number of `ntrials` reduces the fluctuations in the cross section,
        thus, making it more smooth. Smaller number of `ntrials` makes calculations of
        cross section faster.
        """
        return getattr(self, "_glauber_trials", 1000)

    @glauber_trials.setter
    def glauber_trials(self, ntrials):
        self._glauber_trials = ntrials

    def _set_kinematics(self, kin):
        self._projectile_id = {
            lp.pi_plus.pdgid: 1,
            lp.K_plus.pdgid: 3,
            lp.K_S_0.pdgid: 3,
            lp.K_L_0.pdgid: 3,
        }.get(
            abs(kin.p1), 2
        )  # 2 is correct for nucleons and nuclei
        self._lib.xxaini(kin.elab, self._projectile_id, kin.p1.A or 1, kin.p2.A)

    def _generate(self):
        self._lib.psconf()
        # Convert QGSJET to HEPEVT
        self._lib.chepevt()
        return True


class QGSJet2Run(QGSJetRun):
    """Implements all abstract attributes of MCRun for the
    QGSJET-II-xx series of event generators."""

    _event_class = QGSJET2Event

    def _tabulated_cross_section(self, kin=None):
        """Tabulated inelastic or production cross section for
        QGSJET-II-xx event generators."""
        kin = self.kinematics if kin is None else kin
        projectile_id = {
            lp.pi_plus.pdgid: 1,
            lp.K_plus.pdgid: 3,
            lp.K_S_0.pdgid: 3,
            lp.K_L_0.pdgid: 3,
        }.get(
            abs(kin.p1), 2
        )  # 2 is correct for nuclei
        prod = self._lib.qgsect(
            kin.elab,
            projectile_id,
            kin.p1.A or 1,
            kin.p2.A,
        )
        return CrossSectionData(prod=prod)

    def _cross_section(self, kin=None, max_info=False):
        """
        Calculate cross section for QGSJET-II-03/4 event generators.

        Parameters
        ----------
        kin : EventKinematics, optional
        max_info : bool, optional
            Flag indicating whether to calculate the cross section for
            AA collisions using MC. Returns the cross sections for
            specific processes. Default is False and only inelastic or
            production cross section is calculated.

        Returns
        -------
        CrossSectionData
            An object containing the calculated cross section data.

        Notes
        -----
        This method calculates the cross section for the QGSJET-II-0X
        event generator. The cross section is calculated based on the
        provided kinematics and collision type. If `max_info is True`
        and the projectile particle is a nucleus (A > 1), the cross
        section is calculated using the Glauber model. Otherwise, the
        cross section is calculated using the tabulated cross section.
        For h-A collisions, the cross section is calculated using the
        tabulated cross section due to an issue in calculating it from
        QGSJET directly.
        """
        kin = self.kinematics if kin is None else kin

        if max_info and (kin.p1.is_nucleus and kin.p1.A > 1):
            gtot, gprod, _, gdd, gqel, gcoh = self._lib.qgcrossc(self.glauber_trials)
            return CrossSectionData(
                total=gtot,
                prod=gprod,
                diffractive_xb=gdd,
                quasielastic=gqel,
                coherent=gcoh,
            )
        elif (kin.p1.is_nucleus and kin.p1.A > 1) or (
            kin.p2.is_nucleus and kin.p2.A > 1
        ):
            # Note: this is a workaround until the calculated h-A cross
            # sections work correctly
            return self._tabulated_cross_section(kin)

        gtot, gin, gel, gdp, gdt, bel = self._lib.cqgshh_ha_cs()
        return CrossSectionData(
            total=gtot,
            inelastic=gin,
            elastic=gel,
            diffractive_xb=gdp,
            diffractive_ax=gdt,
            b_elastic=bel,
        )

    @property
    def glauber_trials(self):
        """Number of trials for Glauber model integration

        Default is 1000 (set at model initialisation).
        Larger number of `ntrials` reduces the fluctuations in the cross section,
        thus, making it more smooth. Smaller number of `ntrials` makes calculations of
        cross section faster.
        """
        return getattr(self, "_glauber_trials", 1000)

    @glauber_trials.setter
    def glauber_trials(self, ntrials):
        self._glauber_trials = ntrials

    def _set_kinematics(self, kin):
        self._projectile_id = {
            lp.pi_plus.pdgid: 1,
            lp.pi_minus.pdgid: -1,
            lp.proton.pdgid: 2,
            lp.antiproton.pdgid: -2,
            lp.neutron.pdgid: 3,
            lp.antineutron.pdgid: -3,
            lp.K_plus.pdgid: 4,
            lp.K_minus.pdgid: -4,
            lp.K_S_0.pdgid: -5,
            lp.K_L_0.pdgid: 5,
        }.get(
            kin.p1, 2
        )  # 2 is correct for nuclei
        self._lib.qgini(kin.elab, self._projectile_id, kin.p1.A or 1, kin.p2.A)

    def _generate(self):
        self._lib.qgconf()
        # Convert QGSJET to HEPEVT
        self._lib.chepevt()
        return True


class QGSJet01d(QGSJet1Run):
    _version = "01d"
    _library_name = "_qgs01"


class QGSJetII03(QGSJet2Run):
    _version = "II-03"
    _library_name = "_qgsII03"


class QGSJetII04(QGSJet2Run):
    _version = "II-04"
    _library_name = "_qgsII04"
