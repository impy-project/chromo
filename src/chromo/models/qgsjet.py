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

    # needed to skip set_final_state_particles() in MCRun
    _set_final_state_particles_called = True

    def __init__(self, evt_kin, *, seed=None):
        import chromo

        super().__init__(seed)

        # logging
        lun = 6  # stdout
        datdir = _cached_data_dir(self._data_url)
        self._lib.cqgsini(datdir, lun, chromo.debug_level)

        self.kinematics = evt_kin
        self._activate_decay_handler(on=True)
        self._set_final_state_particles()

    def _set_stable(self, pdgid, stable):
        # use Pythia8 instance to decay particles which QGSJet does not decay
        pass


class QGSJet1Run(QGSJetRun):
    """Implements all abstract attributes of MCRun for the
    QGSJET-01c legacy event generators."""

    _event_class = QGSJET1Event

    def _cross_section(self, kin=None):
        # Interpolation routine for QGSJET01D cross sections from CORSIKA.
        from scipy.interpolate import UnivariateSpline

        kin = self.kinematics if kin is None else kin

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

        inel = np.exp(spl(np.log(kin.elab)))
        return CrossSectionData(inelastic=inel)

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

    def _cross_section(self, kin=None):
        kin = self.kinematics if kin is None else kin

        inel = self._lib.qgsect(
            kin.elab,
            self._projectile_id,
            kin.p1.A,
            kin.p2.A,
        )
        return CrossSectionData(inelastic=inel)

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
