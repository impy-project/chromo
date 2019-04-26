'''
Created on 17.03.2014

@author: afedynitch
'''

# Some notes on UrQMD. The model works quite well now. It is slower than anything else,
# except EPOS maybe. It describes quite okaish the fixed target results and some LHC results.
# What is a strange feature is that decays of pions, kaons or neutrinos are not supported
# in the model and can not be disabled by flags.

# The current settings are taken from CORSIKA and they are optimized for speed aparently.
# The license of UrQMD is quite restrictive, they won't probably permit distributing it.



import numpy as np
from impy.common import MCRun, MCEvent, impy_config
from impy.util import standard_particles, info


class UrQMDEvent(MCEvent):
    """Wrapper class around EPOS particle stack."""

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
        raise Exception('Children info not available.')
        # if self._is_filtered:
        #     raise Exception(
        #         'Parent indices do not point to the' +
        #         ' proper particles if any slicing/filtering is applied.')
        # return self.lib.hepevt.jdahep

    @property
    def charge(self):
        return self.lib.uqchg.ichg[self.selection]

    # Nuclear collision parameters
    @property
    def impact_parameter(self):
        """Returns impact parameter for nuclear collisions."""
        return self.lib.rsys.bimp


class UrQMDRun(MCRun):
    """Implements all abstract attributes of MCRun for the 
    UrQMD series of event generators."""

    def __init__(self, *args, **kwargs):
        from particletools.tables import UrQMDParticleTable
        self.stab = UrQMDParticleTable()
        MCRun.__init__(self, *args, **kwargs)

    def sigma_inel(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        return self.lib.ptsigtot()

    def set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""
        k = event_kinematics
        self._curr_event_kin = k

        if k.A1 == 1:
            # Special projectile
            self.lib.inputs.prspflg = 1
            self.lib.sys.ap = 1
            self.lib.inputs.spityp[0] = self.stab.pdg2modid[k.p1pdg][0]
            self.lib.inputs.spiso3[0] = self.stab.pdg2modid[k.p1pdg][1]
        else:
            self.lib.inputs.prspflg = 0
            self.lib.sys.ap = k.A1
            self.lib.sys.zp = k.Z1
        
        if k.A2 == 1:
            # Special projectile
            self.lib.inputs.trspflg = 1
            self.lib.sys.at = 1
            self.lib.inputs.spityp[1] = self.stab.pdg2modid[k.p2pdg][0]
            self.lib.inputs.spiso3[1] = self.stab.pdg2modid[k.p2pdg][1]
        else:
            self.lib.inputs.trspflg = 0
            self.lib.sys.at = k.A2
            self.lib.sys.zt = k.Z2

        # Set impact parameter (to be revisited)
        self.lib.rsys.bdist = (self._nucrad(
            self.lib.sys.ap,
            self.lib.options.ctoption[23]) + self._nucrad(
                self.lib.sys.at, self.lib.options.ctoption[23]) +
                                2 * self.lib.options.ctparam[29])

        # Set ebeam, eos = 0, nevents
        self.lib.input2.pbeam = k.plab
        self.lib.inputs.srtflag = 2
        # Unclear what this thing does but should be required for very low energy
        # if k.plab > 4.9:
        #     self.lib.inputs.eos = 0
        # else:
        #     self.lib.inputs.eos = 1
        #     self.lib.options.ctoption[23]=0

        info(5, 'Setting event kinematics')

    def _nucrad(self, AA, ctopt):
        A = abs(AA)
        nucrad = 0
        rho0 = 0.16
        # root mean square radius of nucleus of mass A
        # r_0 corresponding to rho0
        if ctopt >= 1:
            # root mean square radius of nucleus of mass A (Mayer-Kuckuck)
            nucrad = 1.128 * A**(1. / 3.) - 0.89 * A**(-(1. / 3.))
        else:
            r_0 = (0.75 / np.pi / rho0)**(1. / 3.)
            # subtract gaussian tails, for distributing centroids correctly
            nucrad = r_0 * (0.5 * (A + (A**(1. / 3.) - 1.)**3.))**(1. / 3.)
        return nucrad

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

        self._abort_if_already_initialized()

        if seed == 'random':
            seed = randint(1000000, 10000000)
        else:
            seed = int(seed)
        info(5, 'Using seed:', seed)

        self.lib.init_rmmard(seed)

        self.attach_log()
        info(1, 'First initialization')
        self.lib.urqini(self._lun, impy_config['urqmd']['verbosity'])

        # Set default stable
        self._define_default_fs_particles()
        self.set_event_kinematics(event_kinematics)

        self.lib.sys.eos = 0
        self.lib.inputs.nevents = 1
        self.lib.rsys.bmin = 0
        self.lib.options.ctoption[4] = 1
        # Disable elastic collision
        self.lib.options.ctoption[6] = 1
        

        # Change CTParams and/or CTOptions if needed
        if 'CTParams' in impy_config['urqmd']:
            for ctp in impy_config['urqmd']['CTParams']:
                self.lib.options.ctparams[ctp[0]] = ctp[1]
                info(5, 'CTParams[{}] changed to {}'.format(ctp[0], ctp[1]))
        if 'CTOptions' in impy_config['urqmd']:
            for cto in impy_config['urqmd']['CTOptions']:
                self.lib.options.ctoptions[cto[0]] = cto[1]
                info(5, 'CTOptions[{}] changed to {}'.format(cto[0], cto[1]))

        # Just copied CORSIKA here
        caltim = impy_config['urqmd']['caltim']
        outtim = impy_config['urqmd']['outtim']
        self.lib.pots.dtimestep = outtim
        self._nsteps = int(0.01 + caltim / outtim)
        self._outsteps = 1.
        self.lib.sys.nsteps = self._nsteps
        self.lib.inputs.outsteps = self._outsteps

    def set_stable(self, pdgid, stable=True):
        stable_ids = self.lib.stables.stabvec
        nstab = self.lib.stables.nstable
        try:
            uid = self.stab.pdg2modid[pdgid][0]
        except KeyError:
            info(5, 'Particle {0} unknown to UrQMD'.format(pdgid))
            return
        if stable:
            if uid not in stable_ids[:nstab]:
                self.lib.stables.stabvec[nstab] = uid
                self.lib.stables.nstable += 1
                info(5, 'defining', pdgid, 'as stable particle')
        else:
            if uid in stable_ids[:nstab]:
                self.lib.stables.stabvec[:nstab - 1] = np.array(
                    [pid for pid in stable_ids[:nstab] if pid != uid])
                self.lib.stables.nstable -= 1
                info(5, pdgid, 'allowed to decay')

    def generate_event(self):
        # If new event, initialize projectile and target
        if self.lib.options.ctoption[39] == 0:
            if self.lib.inputs.prspflg == 0:
                self.lib.cascinit(self.lib.sys.zp, self.lib.sys.ap, 1)
            if self.lib.inputs.trspflg == 0:
                self.lib.cascinit(self.lib.sys.zt, self.lib.sys.at, 2)

        self.lib.urqmd(0)
        # Convert URQMD event to HEPEVT
        self.lib.chepevt()
        return 0
