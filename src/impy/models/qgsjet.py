import numpy as np
from impy.common import MCRun, MCEvent
from impy import impy_config
from impy.util import info, _cached_data_dir


class QGSJETEvent(MCEvent):
    """Wrapper class around QGSJet HEPEVT converter."""

    def _charge_init(self, npart):
        return self._lib.qgchg.ichg[:npart]

    @property
    def impact_parameter(self):
        """Returns impact parameter for nuclear collisions."""
        return self._lib.qgarr7.b

    @property
    def n_wounded_A(self):
        """Number of wounded nucleons side A"""
        return self._lib.qgarr55.nwp

    @property
    def n_wounded_B(self):
        """Number of wounded nucleons (target) side B"""
        return self._lib.qgarr55.nwt

    @property
    def n_wounded(self):
        """Number of total wounded nucleons"""
        return self._lib.qgarr55.nwp + self._lib.qgarr55.nwt

    @property
    def n_spectator_A(self):
        """Number of spectator nucleons side A"""
        return self._lib.qgarr56.nspec

    @property
    def n_spectator_B(self):
        """Number of spectator nucleons (target) side B"""
        return self._lib.qgarr56.nspect

    @property
    def diffr_type(self):
        """Type of diffration"""
        return self._lib.jdiff.jdiff


#: Projectiles for QGSJET01 and cross sections
_qgsjet01_projectiles = {211: 1, 111: 1, 2212: 2, 2112: 2, 321: 3, 130: 3, 310: 3}

#: Specific projectile particle indices for QGSII
_qgsjetII_projectiles = {211: 1, 2212: 2, 2112: 3, 321: 4, 130: 5, 310: -5}

#: Used for cross section routines
_qgsjet_hadron_classes = _qgsjet01_projectiles


class QGSJetRun(MCRun):
    _name = "QGSJet"
    _event_class = QGSJETEvent
    _output_frame = "laboratory"
    # needed to skip set_final_state_particles()
    _set_final_state_particles_called = True

    def _set_stable(self, pdgid, stable):
        import warnings

        warnings.warn(
            f"stable particles cannot be changed in {self.pyname}", RuntimeWarning
        )


class QGSJetIIRun(QGSJetRun):
    """Implements all abstract attributes of MCRun for the
    QGSJET-II-xx series of event generators."""

    def __init__(self, event_kinematics, seed=None, logfname=None):

        super().__init__(seed, logfname)

        info(5, "Initializing QGSJET-II")

        datdir = _cached_data_dir(impy_config["qgsjet"]["datdir"])
        self._lib.cqgsini(
            self._seed, datdir, self._lun, impy_config["qgsjet"]["debug_level"]
        )

        self.event_kinematics = event_kinematics

    def _sigma_inel(self, evt_kin):
        return self._lib.qgsect(
            evt_kin.elab,
            _qgsjet_hadron_classes[abs(evt_kin.p1pdg)],
            evt_kin.A1,
            evt_kin.A2,
        )

    def _set_event_kinematics(self, k):
        info(5, "Setting event kinematics")

        if abs(k.p1pdg) in _qgsjetII_projectiles:
            self._qgsproj = np.sign(k.p1pdg) * _qgsjetII_projectiles[abs(k.p1pdg)]
        elif k.p1pdg == 111:
            # For pi0 projectiles alternate between pi+ and pi-
            self._qgsproj = int(1 - 2 * round(np.random.rand()))
        else:
            raise Exception(
                "Projectile {0} not supported by QGSJET-II.".format(k.p1pdg)
            )

        self._lib.qgini(k.elab, self._qgsproj, k.A1, k.A2)

    def _attach_log(self, fname=None):
        """Routes the output to a file or the stdout."""
        fname = impy_config["output_log"] if fname is None else fname
        if fname == "stdout":
            lun = 6
            info(5, "Output is routed to stdout.")
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, "Output is routed to", fname, "via LUN", lun)
        self._lun = lun

    def _generate_event(self):
        self._lib.qgconf()
        # Convert QGSJET to HEPEVT
        self._lib.chepevt()
        return False


class QGSJet01Run(QGSJetRun):
    """Implements all abstract attributes of MCRun for the
    QGSJET-01c legacy event generators."""

    def __init__(self, event_kinematics, seed="random", logfname=None):
        super().__init__(seed, logfname)

        info(5, "Initializing QGSJET01d")

        datdir = _cached_data_dir(impy_config["qgsjet"]["datdir"])
        self._lib.cqgsini(
            self._seed, datdir, self._lun, impy_config["qgsjet"]["debug_level"]
        )

        self.event_kinematics = event_kinematics

    def _sigma_inel(self, evt_kin):
        # Interpolation routine for QGSJET01D cross sections from CORSIKA.
        from scipy.interpolate import UnivariateSpline

        A_target = evt_kin.A2
        # Projectile ID-1 to access fortran indices directly
        icz = self._qgsproj - 1
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

        return np.exp(spl(np.log(evt_kin.elab)))

    def _set_event_kinematics(self, k):
        info(5, "Setting event kinematics")

        if abs(k.p1pdg) in [2212, 2112]:
            self._qgsproj = 2
        elif abs(k.p1pdg) in [211, 111]:
            self._qgsproj = 1
        elif abs(k.p1pdg) in [321, 130, 310]:
            self._qgsproj = 3
        else:
            raise Exception("QGSJET only supports p, pi+- and K+- as projectile.")
        self._lib.xxaini(k.elab, self._qgsproj, k.A1, k.A2)

    def _attach_log(self, fname=None):
        """Routes the output to a file or the stdout."""
        fname = impy_config["output_log"] if fname is None else fname
        if fname == "stdout":
            lun = 6
            info(5, "Output is routed to stdout.")
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, "Output is routed to", fname, "via LUN", lun)

        self._lun = lun

    def _generate_event(self):
        self._lib.psconf()
        # Convert QGSJET to HEPEVT
        self._lib.chepevt()
        return False


class QGSJet01d(QGSJet01Run):
    _version = "01d"
    _library_name = "_qgs01"


class QGSJetII03(QGSJetIIRun):
    _version = "II-03"
    _library_name = "_qgsII03"


class QGSJetII04(QGSJetIIRun):
    _version = "II-04"
    _library_name = "_qgsII04"
