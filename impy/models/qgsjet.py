'''
Created on 17.03.2014

@author: afedynitch
'''
import numpy as np
from impy.common import MCRun, MCEvent, impy_config, root_dir
from impy.util import standard_particles, info

class QGSJETEvent(MCEvent):
    """Wrapper class around QGSJet HEPEVT converter."""

    def __init__(self, lib, event_kinematics, event_frame):
        # HEPEVT (style) common block
        evt = lib.hepevt

        # Save selector for implementation of on-demand properties
        px, py, pz, en, m = evt.phep
        vx, vy, vz, vt = evt.vhep

        MCEvent.__init__(
            self,
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
        if self._is_filtered:
            raise Exception(
                'Parent indices do not point to the' +
                ' proper particles if any slicing/filtering is applied.')
        return self.lib.hepevt.jmohep

    @property
    def children(self):
        if self._is_filtered:
            raise Exception(
                'Parent indices do not point to the' +
                ' proper particles if any slicing/filtering is applied.')
        return self.lib.hepevt.jdahep

    @property
    def charge(self):
        return self.lib.qgchg.ichg[self.selection]

class QGSJetIIRun(MCRun):
    """Implements all abstract attributes of MCRun for the 
       QGSJET-II-xx series of event generators."""

    def __init__(self, *args, **kwargs):
        from particletools.tables import QGSJetParticleTable
        self.stab = QGSJetParticleTable()
        MCRun.__init__(self, *args, **kwargs)

    def sigma_inel(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        k = self._curr_event_kin
        return self.lib.qgsect(self._curr_event_kin.elab, self._qgsproj, k.A1,
                               k.A2)

    def hadron_air_cs(self):
        """Hadron-air production cross sections according to current
        event setup (energy, projectile)."""
        # Mass composition of air (Nitrogen, Oxygen, Argon)
        frac_air = [(0.78479, 14), (0.21052, 16), (0.00469, 40)]
        return np.sum([
            f * self.lib.qgsect(self._curr_event_kin.elab, self._qgsproj, 1,
                                iat) for f, iat in frac_air
        ], axis=0)

    def set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""

        info(5, 'Setting event kinematics')
        k = event_kinematics
        self._curr_event_kin = k
        self._qgsproj = abs(self.stab.pdg2modid[k.p1pdg])
        # if k.p1pdg not in self.impy_config['qgsjet']

        if not 0 < self._qgsproj < 4:
            raise Exception(
                'QGSJET only supports p, pi+- and K+- as projectile.')
        self.lib.qgini(k.elab, self._qgsproj, k.A1, k.A2)

    def attach_log(self):
        """Routes the output to a file or the stdout."""
        fname = impy_config['output_log']
        if fname == 'stdout':
            lun = 6
            info(5, 'Output is routed to stdout.')
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, 'Output is routed to', fname, 'via LUN', lun)

        self._lun = lun

    def init_generator(self, event_kinematics, seed='random'):
        from random import randint
        from os import path

        self._abort_if_already_initialized()

        if seed == 'random':
            seed = randint(1000000, 10000000)
        else:
            seed = int(seed)
        info(5, 'Using seed:', seed)

        info(5, 'Initializing QGSJET-II')
        datdir = path.join(root_dir, impy_config['qgsjet']['datdir'])
        self.attach_log()
        self.lib.cqgsini(seed, datdir, self._lun)

        # Set default stable
        info(10, 'All particles stable in QGSJET-II')
        self.set_event_kinematics(event_kinematics)
    
    def set_stable(self, stable):
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
        self.stab = QGSJetParticleTable()
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
        k = self._curr_event_kin
        info(2,'Sigma inel not implemented for QGSJet01c')
        return 0.

    def hadron_air_cs(self):
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

        spl = UnivariateSpline(
            np.log(qgsgrid),
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
        self._qgsproj = abs(self.stab.pdg2modid[k.p1pdg])

        if not 0 < self._qgsproj < 4:
            raise Exception(
                'QGSJET only supports p, pi+- and K+- as projectile.')
        self.lib.xxaini(k.elab, self._qgsproj, k.A1, k.A2)

    def attach_log(self):
        """Routes the output to a file or the stdout."""
        fname = impy_config['output_log']
        if fname == 'stdout':
            lun = 6
            info(5, 'Output is routed to stdout.')
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, 'Output is routed to', fname, 'via LUN', lun)

        self._lun = lun

    def init_generator(self, event_kinematics, seed='random'):
        from random import randint
        from os import path

        self._abort_if_already_initialized()

        if seed == 'random':
            seed = randint(1000000, 10000000)
        else:
            seed = int(seed)
        info(5, 'Using seed:', seed)

        info(5, 'Initializing QGSJET01c')
        self.attach_log()
        datdir = path.join(root_dir, impy_config['qgsjet']['datdir'])
        self.lib.cqgsini(seed, datdir, self._lun)

        # Set default stable
        info(10, 'All particles stable in QGSJET-01')
        self.set_event_kinematics(event_kinematics)

    def set_stable(self, stable):
        info(10, "All particles stable in QGSJet.")

    def generate_event(self):
        self.lib.psconf()
        # Convert QGSJET to HEPEVT
        self.lib.chepevt()
        return False