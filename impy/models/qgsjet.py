'''
Created on 17.03.2014

@author: afedynitch
'''
import numpy as np
from impy.common import MCRun, MCEvent
from impy import impy_config, base_path
from impy.util import standard_particles, info


class QGSJETEvent(MCEvent):
    """Wrapper class around QGSJet HEPEVT converter."""
    def __init__(self, lib, event_kinematics, event_frame):
        # HEPEVT (style) common block
        evt = lib.hepevt

        # Save selector for implementation of on-demand properties
        px, py, pz, en, m = evt.phep
        vx, vy, vz, vt = evt.vhep

        MCEvent.__init__(self,
                         lib=lib,
                         event_kinematics=event_kinematics,
                         event_frame=event_frame,
                         nevent=evt.nevhep,
                         npart=evt.nhep,
                         p_ids=evt.idhep,
                         status=evt.isthep,
                         px=px,
                         py=py,
                         pz=pz,
                         en=en,
                         m=m,
                         vx=vx,
                         vy=vy,
                         vz=vz,
                         vt=vt,
                         pem_arr=evt.phep,
                         vt_arr=evt.vhep)

    def filter_final_state(self):
        self.selection = np.where(self.status == 1)
        self._apply_slicing()

    def filter_final_state_charged(self):
        self.selection = np.where((self.status == 1) & (self.charge != 0))
        self._apply_slicing()

    @property
    def parents(self):
        MCEvent.parents(self)
        return self.lib.hepevt.jmohep

    @property
    def children(self):
        MCEvent.children(self)
        return self.lib.hepevt.jdahep

    @property
    def charge(self):
        return self.lib.qgchg.ichg[self.selection]

    @property
    def impact_parameter(self):
        """Returns impact parameter for nuclear collisions."""
        return self.lib.qgarr7.b

    @property
    def n_wounded_A(self):
        """Number of wounded nucleons side A"""
        return self.lib.qgarr55.nwp

    @property
    def n_wounded_B(self):
        """Number of wounded nucleons (target) side B"""
        return self.lib.qgarr55.nwt

    @property
    def n_wounded(self):
        """Number of total wounded nucleons"""
        return self.lib.qgarr55.nwp + self.lib.qgarr55.nwt

    @property
    def n_spectator_A(self):
        """Number of spectator nucleons side A"""
        return self.lib.qgarr56.nspec

    @property
    def n_spectator_B(self):
        """Number of spectator nucleons (target) side B"""
        return self.lib.qgarr56.nspect

    @property
    def diffr_type(self):
        """Type of diffration"""
        return self.lib.jdiff.jdiff


#: Projectiles for QGSJET01 and cross sections
_qgsjet01_projectiles = {
    211: 1,
    111: 1,
    2212: 2,
    2112: 2,
    321: 3,
    130: 3,
    310: 3
}

#: Specific projectile particle indices for QGSII
_qgsjetII_projectiles = {211: 1, 2212: 2, 2112: 3, 321: 4, 130: 5, 310: -5}

#: Used for cross section routines
_qgsjet_hadron_classes = _qgsjet01_projectiles


class QGSJetIIRun(MCRun):
    """Implements all abstract attributes of MCRun for the 
       QGSJET-II-xx series of event generators."""
    def __init__(self, *args, **kwargs):
        from particletools.tables import QGSJetParticleTable
        MCRun.__init__(self, *args, **kwargs)

    def sigma_inel(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        k = self._curr_event_kin
        return self.lib.qgsect(self._curr_event_kin.elab,
                               _qgsjet_hadron_classes[abs(k.p1pdg)], k.A1,
                               k.A2)

    def set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target for next event."""

        info(5, 'Setting event kinematics')
        k = event_kinematics
        self._curr_event_kin = k
        if abs(k.p1pdg) in _qgsjetII_projectiles:
            self._qgsproj = np.sign(k.p1pdg) * _qgsjetII_projectiles[abs(
                k.p1pdg)]
        elif k.p1pdg == 111:
            # For pi0 projectiles alternate between pi+ and pi-
            self._qgsproj = int(1 - 2 * round(np.random.rand()))
        else:
            raise Exception(
                'Projectile {0} not supported by QGSJET-II.'.format(k.p1pdg))

        self.lib.qgini(k.elab, self._qgsproj, k.A1, k.A2)

    def attach_log(self, fname=None):
        """Routes the output to a file or the stdout."""
        fname = impy_config['output_log'] if fname is None else fname
        if fname == 'stdout':
            lun = 6
            info(5, 'Output is routed to stdout.')
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, 'Output is routed to', fname, 'via LUN', lun)
        self._lun = lun

    def init_generator(self, event_kinematics, seed='random', logfname=None):
        from random import randint
        from os import path

        self._abort_if_already_initialized()

        if seed == 'random':
            seed = randint(1000000, 10000000)
        else:
            seed = int(seed)
        info(5, 'Using seed:', seed)

        info(5, 'Initializing QGSJET-II')
        datdir = path.join(base_path, impy_config['qgsjet']['datdir'])
        self.attach_log(fname=logfname)
        self.lib.cqgsini(seed, datdir, self._lun,
                         impy_config['qgsjet']['debug_level'])

        # Set default stable
        info(10, 'All particles stable in QGSJET-II')
        self.set_event_kinematics(event_kinematics)

    def set_stable(self, pdgid, stable=True):
        info(10, "All particles stable in QGSJet-II.")

    def generate_event(self):
        self.lib.qgconf()
        # Convert QGSJET to HEPEVT
        self.lib.chepevt()
        return False


class QGSJet01Run(MCRun):
    """Implements all abstract attributes of MCRun for the 
       QGSJET-01c legacy event generators."""
    def __init__(self, *args, **kwargs):
        from particletools.tables import QGSJetParticleTable
        MCRun.__init__(self, *args, **kwargs)

    # def sigma_inel(self):
    #     """Inelastic cross section according to current
    #     event setup (energy, projectile, target)"""
    #     k = self._curr_event_kin
    #     return self.lib.qgsect(self._curr_event_kin.elab, self._qgsproj, k.A1,
    #                            k.A2)

    def sigma_inel(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        # k = self._curr_event_kin
        info(2, 'Sigma inel not implemented for QGSJet01c')
        return 0.

    def sigma_inel_air(self):
        """Hadron-air production cross sections according to current
        event setup (energy, projectile)."""
        # Mass composition of air (Nitrogen, Oxygen, Argon)
        from scipy.interpolate import UnivariateSpline

        # Projectile ID-1 to access fortran indices directly
        icz = self._qgsproj - 1

        qgsgrid = 10**np.arange(1, 11)
        cross_section = np.zeros(10)
        frac_air = [(0.78479, 14), (0.21052, 16), (0.00469, 40)]
        wa = np.zeros(3)
        for frac, iat in frac_air:
            for je in range(10):
                sectn = 0.

                ya = iat

                ya = np.log(ya) / 1.38629 + 1.
                ja = min(int(ya), 2)
                wa[1] = (ya - ja)
                wa[2] = wa[1] * (wa[1] - 1) * .5
                wa[0] = 1. - wa[1] + wa[2]
                wa[1] = wa[1] - 2. * wa[2]
                for m in range(3):
                    sectn += self.lib.xsect.gsect[je, icz, ja + m - 1] * wa[m]

                cross_section[je] += frac * np.exp(sectn)

        spl = UnivariateSpline(np.log(qgsgrid),
                               np.log(cross_section),
                               ext='extrapolate',
                               s=0,
                               k=1)
        return np.exp(spl(np.log(self._curr_event_kin.elab)))

    def set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""

        info(5, 'Setting event kinematics')
        k = event_kinematics
        self._curr_event_kin = k

        if abs(k.p1pdg) in [2212, 2112]:
            self._qgsproj = 2
        elif abs(k.p1pdg) in [211, 111]:
            self._qgsproj = 1
        elif abs(k.p1pdg) in [321, 130, 310]:
            self._qgsproj = 3
        else:
            raise Exception(
                'QGSJET only supports p, pi+- and K+- as projectile.')
        self.lib.xxaini(k.elab, self._qgsproj, k.A1, k.A2)

    def attach_log(self, fname=None):
        """Routes the output to a file or the stdout."""
        fname = impy_config['output_log'] if fname is None else fname
        if fname == 'stdout':
            lun = 6
            info(5, 'Output is routed to stdout.')
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, 'Output is routed to', fname, 'via LUN', lun)

        self._lun = lun

    def init_generator(self, event_kinematics, seed='random', logfname=None):
        from random import randint
        from os import path

        self._abort_if_already_initialized()

        if seed == 'random':
            seed = randint(1000000, 10000000)
        else:
            seed = int(seed)
        info(5, 'Using seed:', seed)

        info(5, 'Initializing QGSJET01c')
        self.attach_log(fname=logfname)
        datdir = path.join(base_path, impy_config['qgsjet']['datdir'])
        self.lib.cqgsini(seed, datdir, self._lun,
                         impy_config['qgsjet']['debug_level'])

        # Set default stable
        info(10, 'All particles stable in QGSJET-01')
        self.set_event_kinematics(event_kinematics)

    def set_stable(self, pdgid, stable=True):
        info(10, "All particles stable in QGSJet.")

    def generate_event(self):
        self.lib.psconf()
        # Convert QGSJET to HEPEVT
        self.lib.chepevt()
        return False