'''
Created on 17.03.2014

@author: afedynitch
'''

import numpy as np
from impy.common import MCRun, MCEvent, impy_config 
from impy.util import standard_particles, info
from particletools.tables import UrQMDParticleTable

class UrQMDMCEvent(MCEvent):
    def __init__(self, lib, event_config):
        # Number of entries on stack
        npart = lib.sys.npart
        # Particle IDs for UrQMD
        list_particle_ids = np.column_stack((lib.isys.ityp, lib.isys.iso3))[:npart]
        # UrQMD to PDG conversion (using table)
        stab = UrQMDParticleTable()
        to_pdg = np.vectorize(lambda x: stab.modid2pdg(tuple(x)))

        # Filter stack for stable and/or charged particles if selected
        stable = ([i for i, p in enumerate(list_particle_ids) 
                   if p[0] in lib.stables.stabvec])
        if impy_config["event_scope"] == 'charged':
            sel = stable[lib.isys.charge[stable] != 0]
        elif impy_config["event_scope"] == 'stable':
            sel = stable
        else:
            raise Exception("not implemented, yet")

        # Save selector for implementation of on-demand properties
        self.sel = sel

        MCEvent.__init__(
            self,
            event_config=event_config,
            lib=lib,
            px=lib.coor.px[sel],
            py=lib.coor.py[sel],
            pz=lib.coor.pz[sel],
            en=lib.coor.p0[sel],
            p_ids=to_pdg(list_particle_ids[sel]),
            npart=npart)

    @property
    def charge(self):
        self.charge = self.lib.isys.charge[self.sel]

    @property
    def mass(self):
        return self.lib.coor.fmass

    # Nuclear collision parameters

    @property
    def impact_parameter(self):
        """Impact parameter for nuclear collisions."""
        return self.lib.rsys.bimp


class UrQMDMCRun(MCRun):
    """Implements all abstract attributes of MCRun for 
    UrQMD 3.4."""

    def __init__(self, libref, event_class=None, **kwargs):
        if event_class is None:
            self.event_class = UrQMDMCEvent
        else:
            self.event_class = event_class

        self._frame = 'center-of-mass'

        MCRun.__init__(self, libref, **kwargs)
    
    @property
    def frame(self):
        return self._frame

    @property
    def name(self):
        """Event generator name"""
        return "UrQMD"

    @property
    def version(self):
        """Event generator version"""
        return "3.4"

    def sigma_inel(self):
        raise Exception("Inelastic cross section not implemented yet")

    def set_event_kinematics(self, event_kinematics):
        stab = UrQMDParticleTable()
        # Define projectile
        if event_kinematics.A1 == 1:
            # Special projectile
            self.lib.inputs.prspflg = 1
            self.lib.sys.ap = 1
            self.lib.inputs.spityp[0] = stab.pdg2modid[event_kinematics.p1pdg][0]
            self.lib.inputs.spiso3[0] = stab.pdg2modid[event_kinematics.p1pdg][1]
        else:
            # Nucleus
            self.lib.inputs.prspflg = 0
            self.lib.sys.ap = event_kinematics.A1
            self.lib.sys.zp = event_kinematics.Z1

        # Define target
        if event_kinematics.A2 == 1:
            # Special target
            self.lib.inputs.trspflg = 1
            self.lib.sys.at = 1
            self.lib.inputs.spityp[1] = stab.pdg2modid[event_kinematics.p2pdg][0]
            self.lib.inputs.spiso3[1] = stab.pdg2modid[event_kinematics.p2pdg][1]
        else:
            # Nucleus
            self.lib.inputs.trspflg = 0
            self.lib.sys.at = event_kinematics.A2
            self.lib.sys.zt = event_kinematics.Z2
        # Set ebeam, eos = 0, nevents
        self.lib.rsys.ebeam = event_kinematics.elab # Needs to be fixed!!
        self.lib.sys.eos = 0
        self.lib.inputs.nevents = 1

        # Set impact parameter (to be revisited)
        self.lib.rsys.bmin = 0
        self.lib.options.ctoption[4] = 1
        self.lib.rsys.bdist = nucrad(self.lib.sys.ap, self.lib.options.ctoption[23]) + nucrad(self.lib.sys.at, self.lib.options.ctoption[23]) + 2 * self.lib.options.ctparam[29]
        self._curr_event_kin = event_kinematics

    def attach_log(self, fname):
        """Routes the output to a file or the stdout."""
        info(5, "Not implemented yet!")
        # TO BE IMPLEMENTED!!!
        if fname == 'stdout':
            pass
        else:
            pass

    def init_generator(self, event_kinematics):
        from random import randint

        self._abort_if_already_initialized()

        self.set_event_kinematics(event_kinematics)

        self.attach_log('stdout')

        set_stable(self.lib, 2)


        # if self.def_settings:
            # print(self.class_name + "::init_generator(): Using default settings:", \
                # self.def_settings.__class__.__name__)
            # self.def_settings.enable()

        # REVISIT HERE!!
        # self._define_default_fs_particles()

        # Set pi0 stable
        self.set_stable(111)

        # if impy_config['stable']:
            # for pdgid in impy_config['stable']:
                # self.set_stable(pdgid)

        self.lib.urqini('stdout', 2)

        # Change CTParams and/or CTOptions if needed
        if 'CTParams' in impy_config:
            for ctp in impy_config['CTParams']:
                self.lib.options.ctparams[ctp[0]] = ctp[1]
                info(5, 'CTParams[{}] changed to {}'.format(ctp[0], ctp[1]))
        if 'CTOptions' in impy_config:
            for ctp in impy_config['CTOptions']:
                self.lib.options.ctoptions[ctp[0]] = ctp[1]
                info(5, 'CTOptions[{}] changed to {}'.format(ctp[0], ctp[1]))

    def set_stable(self, pdgid):
        if pdgid == 0:
            return
        info(5, 'defining', pdgid, 'as stable particle')
        stab = UrQMDParticleTable()
        stable_count = np.where(self.lib.stables.stabvec == 0)[0][0]
        self.lib.stables.stabvec[stable_count] = stab.pdg2modid[pdgid][0]
        self.lib.stables.nstable = stable_count + 1

    def generate_event(self, iflbmax):
        # Set calculation time (can be introduced as arguments later)
        # Just copied CORSIKA here
        caltim = 200.0
        outtim = 200.0
        self.lib.pots.dtimestep = outtim
        self.lib.sys.nsteps = int(0.01 + caltim/self.lib.pots.dtimestep)
        self.lib.sys.outsteps = int(0.01 + outtim/self.lib.pots.dtimestep)

        # If new event, initialize projectile and target
        if self.lib.options.ctoption[39] == 0:
            if self.lib.inputs.prspflg == 0:
                self.lib.cascinit(self.lib.sys.zp, self.lib.sys.ap, 1)
            if self.lib.inputs.trspflg == 0:
                self.lib.cascinit(self.lib.sys.zt, self.lib.sys.at, 2)

        self.lib.urqmd(iflbmax)
        return 0

#=========================================================================
# set_stable
#=========================================================================
def set_stable(lib, decay_mode, dbg=True):
    lib.stables.stabvec = np.zeros((len(lib.stables.stabvec),))

    if dbg:
        print("UrQMDCascadeRun::set_stable(): Setting standard" +
              " particles stable.")

    # Count number of stable particles
    stable_count = 0

    #fast-mode particles
    if decay_mode == 0:
        stab = UrQMDParticleTable()
        for pdg_id in standard_particles:
            try:
                lib.stables.stabvec[stable_count] = stab.pdg2modid[pdg_id][0]
                info(5, 'stable,', pdg_id)
                stable_count += 1
            except KeyError:
                pass
        lib.stables.nstable = stable_count
        return
    # STILL TO CODE: decay_mode > 0 cases

    # # keep muons pions, kaons
    # for i in range(4, 5 + 1):
    # idb[i - 1] = -np.abs(idb[i - 1])
    # for i in range(7, 18 + 1):
    # idb[i - 1] = -np.abs(idb[i - 1])
    # # K0 and K0-bar have to remain unstable to form K0S/L

    # if decay_mode <= 1:
    # self.stables.nstable = (idb != 0).sum()
    # return

    # # Decay mode 2 for generation of decay spectra (all conventional with
    # # lifetime >= K0S
    # if dbg:
    # print("UrQMDCascadeRun::set_stable(): Setting conventional " +
    # "Sigma-, Xi0, Xi- and Lambda0 stable (decay mode).")
    # for i in range(36, 39 + 1):
    # idb[i - 1] = -np.abs(idb[i - 1])

    # if decay_mode <= 2:
    # self.stables.nstable = (idb != 0).sum()
    # return

    # # Conventional mesons and baryons
    # # keep eta, eta', rho's, omega, phi, K*
    # if dbg:
    # print("UrQMDCascadeRun::set_stable(): Setting all " +
    # "conventional stable.")
    # # pi0
    # idb[6 - 1] = -np.abs(idb[6 - 1])
    # for i in range(23, 33 + 1):
    # idb[i - 1] = -np.abs(idb[i - 1])

    # # keep SIGMA, XI, LAMBDA
    # for i in range(34, 49 + 1):
    # idb[i - 1] = -np.abs(idb[i - 1])

    # if decay_mode <= 3:
    # self.stables.nstable = (idb != 0).sum()
    # return

    # # Charmed particles (only for version >= 2.2)
    # # keep all charmed
    # if dbg:
    # print("UrQMDCascadeRun::set_stable(): Setting all " +
    # "conventional and charmed stable.")
    # for i in range(59, 61) + range(71, 99 + 1):
    # idb[i - 1] = -np.abs(idb[i - 1])
    # self.stables.nstable = (idb != 0).sum()

def nucrad(AA, ctopt):
    A=abs(AA)
    nucrad = 0
    rho0 = 0.16
    # root mean square radius of nucleus of mass A
    # r_0 corresponding to rho0
    if ctopt >= 1:
        # root mean square radius of nucleus of mass A (Mayer-Kuckuck)
        nucrad = 1.128 * A**(1./3.) - 0.89 * A**(-(1./3.))
    else:
        r_0 = (0.75/np.pi/rho0)**(1./3.)
        # subtract gaussian tails, for distributing centroids correctly
        nucrad = r_0*(0.5*(A + (A**(1./3.)-1.)**3.))**(1./3.)
    return nucrad
