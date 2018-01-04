'''
Created on 17.03.2014

@author: afedynitch
'''

import numpy as np
from impy.common import MCRun, MCEvent, EventKinematics, standard_particles
from ParticleDataTool import UrQMDParticleTable

class UrQMDMCEvent(MCEvent):
    def __init__(self, lib, event_config):
        npart = lib.sys.npart
        self.p_ids = np.zeros((npart,))
        list_particle_ids = np.concatenate((np.atleast_2d(lib.isys.ityp), np.atleast_2d(lib.isys.iso3)), axis = 0)[:npart].transpose()
        stable = ([i for i, p in enumerate(list_particle_ids) if p[0] in lib.stables.stabvec])
        stab = UrQMDParticleTable()
        sel = None
        if event_config['charged_only']:
            sel = stable[lib.isys.charge[stable] != 0]
        else:
            sel = stable

        if 'charge_info' in event_config and event_config['charge_info']:
            self.charge = lib.isys.charge[sel]

        self.p_ids = np.array([ stab.modid2pdg[tuple(l)] 
                            for l in list_particle_ids[stable]])
        self.en = lib.coor.p0[sel]
        self.pz = lib.coor.pz[sel]
        self.pt2 = lib.coor.py[sel]**2 + lib.coor.px[sel]**2

        MCEvent.__init__(self, event_config)


class UrQMDMCRun(MCRun):
    def __init__(self, libref, **kwargs):

        if "event_class" not in kwargs.keys():
            kwargs["event_class"] = UrQMDMCEvent

        MCRun.__init__(self, libref, **kwargs)

        self.generator = "UrQMD"
        self.version = "3.4"
        self.ptab = UrQMDParticleTable()

    def get_sigma_inel(self):
        raise Exception(self.class_name + "::init_generator(): not implemented yet")

    def set_event_kinematics(self):
        # Define projectile
        if self.evkin.A1 == 1:
            # Special projectile
            self.lib.inputs.prspflg = 1
            self.lib.sys.ap = 1
            self.lib.inputs.spityp[0] = self.ptab.pdg2modid[self.evkin.p1pdg][0]
            self.lib.inputs.spiso3[0] = self.ptab.pdg2modid[self.evkin.p1pdg][1]
        else:
            # Nucleus
            self.lib.inputs.prspflg = 0
            self.lib.sys.ap = self.evkin.A1
            self.lib.sys.zp = self.evkin.Z1

        # Define target
        if self.evkin.A2 == 1:
            # Special target
            self.lib.inputs.trspflg = 1
            self.lib.sys.at = 1
            self.lib.inputs.spityp[1] = self.ptab.pdg2modid[self.evkin.p2pdg][0]
            self.lib.inputs.spiso3[1] = self.ptab.pdg2modid[self.evkin.p2pdg][1]
        else:
            # Nucleus
            self.lib.inputs.trspflg = 0
            self.lib.sys.at = self.evkin.A2
            self.lib.sys.zt = self.evkin.Z2
        # Set ebeam, eos = 0, nevents
        self.lib.rsys.ebeam = self.evkin.elab # Needs to be fixed!!
        self.lib.sys.eos = 0
        self.lib.inputs.nevents = 1

        # Set impact parameter (to be revisited)
        self.lib.rsys.bmin = 0
        self.lib.options.ctoption[4] = 1
        self.lib.rsys.bdist = nucrad(self.lib.sys.ap, self.lib.options.ctoption[23]) + nucrad(self.lib.sys.at, self.lib.options.ctoption[23]) + 2 * self.lib.options.ctparam[29]
        self.event_config['event_kinematics'] = self.evkin

    def init_generator(self, config):

        try:
            print self.init
            raise Exception('adasd')
        except AttributeError:
            self.init = True

        self.set_event_kinematics()

        from random import randint

        file_output = ''
        try:
            from MCVD.management import LogManager
            # initialize log manager
            self.log_man = LogManager(config, self)
            self.log_man.create_log(self.sibproj, self.eatarg, self.ecm)
            # Still need to determine where the file is outputted
        except ImportError:
            print self.class_name + "::init_generator(): Running outside of MCVD,", \
                "the log will be printed to STDOUT."
        except AttributeError:
            print self.class_name + "::init_generator(): Logging not supported."
        set_stable(self.lib, 2)
        self.set_event_kinematics()

        if self.def_settings:
            print self.class_name + "::init_generator(): Using default settings:", \
                self.def_settings.__class__.__name__
            self.def_settings.enable()

        # Set pi0 stable
        self.set_stable(111)

        if 'stable' in self.event_config:
            for pdgid in self.event_config['stable']:
                self.set_stable(pdgid)

        self.lib.urqini(file_output, 2)
        # Change CTParams and/or CTOptions if needed
        if hasattr(self, 'CTParams'):
            for ctp in self.CTParams:
                self.lib.options.ctparams[ctp[0]] = ctp[1]
                print self.__class__.__name__ + '::init_generator(): CTParams[{}] changed to {}'.format(ctp[0], ctp[1])
        if hasattr(self, 'CTOptions'):
            for ctp in self.CTOptions:
                self.lib.options.ctoptions[ctp[0]] = ctp[1]
                print self.__class__.__name__ + '::init_generator(): CTOptions[{}] changed to {}'.format(ctp[0], ctp[1])
        print self.__class__.__name__ + '::init_generator(): Done'

    def set_stable(self, pdgid):
        if pdgid == 0:
            return
        print self.class_name + "::set_stable(): defining ", \
            pdgid, "as stable particle, sid =", pdgid
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

class UrQMDCascadeRun():
    def __init__(self,
                 lib_str,
                 label,
                 decay_mode,
                 n_events,
                 fill_subset=False,
                 p_debug=False,
                 nucleon_Ekin=None,
                 CTParams = None,
                 CTOptions = None):
        exec "import " + lib_str + " as urlib"
        self.lib = urlib
        self.label = label
        self.nEvents = n_events
        self.fill_subset = fill_subset
        self.spectrum_hists = []
        self.dbg = p_debug
        self.nucleon_Ekin = nucleon_Ekin
        set_stable(self.lib, decay_mode, self.dbg)
        self.ptab = UrQMDParticleTable()
        self.CTParams = CTParams
        self.CTOptions = CTOptions
        self.init_generator()

    def init_generator(self):
        self.lib.urqini('output_urtest.txt', 2)
        # Change CTParams and/or CTOptions if needed
        if self.CTParams:
            for ctp in self.CTParams:
                self.lib.options.ctparams[ctp[0]] = ctp[1]
                print self.__class__.__name__ + '::init_generator(): CTParams[{}] changed to {}'.format(ctp[0], ctp[1])
        if self.CTOptions:
            for ctp in self.CTOptions:
                self.lib.options.ctoptions[ctp[0]] = ctp[1]
                print self.__class__.__name__ + '::init_generator(): CTOptions[{}] changed to {}'.format(ctp[0], ctp[1])
        print self.__class__.__name__ + '::init_generator(): Done'

    def get_hadron_air_cs(self, E_lab, projectile_sibid):
        raise Exception('Not implemented, yet')

    def start(self, evkin, iflbmax):
        # Define projectile
        if evkin.A1 == 1:
            # Special projectile
            self.lib.inputs.prspflg = 1
            self.lib.sys.ap = 1
            self.lib.inputs.spityp[0] = self.ptab.pdg2modid[evkin.p1pdg][0]
            self.lib.inputs.spiso3[0] = self.ptab.pdg2modid[evkin.p1pdg][1]
        else:
            # Nucleus
            self.lib.inputs.prspflg = 0
            self.lib.sys.ap = evkin.A1
            self.lib.sys.zp = evkin.Z1

        # Define target
        if evkin.A2 == 1:
            # Special target
            self.lib.inputs.trspflg = 1
            self.lib.sys.at = 1
            self.lib.inputs.spityp[1] = self.ptab.pdg2modid[evkin.p2pdg][0]
            self.lib.inputs.spiso3[1] = self.ptab.pdg2modid[evkin.p2pdg][1]
        else:
            # Nucleus
            self.lib.inputs.trspflg = 0
            self.lib.sys.at = evkin.A2
            self.lib.sys.zt = evkin.Z2
        # Set ebeam, eos = 0, nevents
        self.lib.rsys.ebeam = evkin.elab
        self.lib.sys.eos = 0
        self.lib.inputs.nevents = 1

        # Set impact parameter (to be revisited)
        self.lib.rsys.bmin = 0
        self.lib.options.ctoption[4] = 1
        self.lib.rsys.bdist = nucrad(self.lib.sys.ap, self.lib.options.ctoption[23]) + nucrad(self.lib.sys.at, self.lib.options.ctoption[23]) + 2 * self.lib.options.ctparam[29]

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

        hist_d = {}
        for hist in self.spectrum_hists:
            hist_d[hist.particle_id] = hist
        ngenerated = self.nEvents

        for i in xrange(self.nEvents):
            self.lib.urqmd(iflbmax)
            if not (i % 10000) and i and self.dbg:
                print i, "events generated."

            event = UrQMDCascadeEvent(self.lib, self.ptab)

            unique_pids = np.unique(event.p_ids)
            
            if 0 in unique_pids:
                ngenerated = ngenerated - 1
                continue
            if not self.fill_subset:
                [hist_d[pid].fill_event(event) for pid in unique_pids]
            else:
                for pid in unique_pids:
                    if pid in hist_d.keys():
                        hist_d[pid].fill_event(event)
        # Correct for selective filling of histograms
        for hist in self.spectrum_hists:
            hist.n_events_filled = ngenerated


class UrQMDCascadeEvent():
    def __init__(self, lib, ptab, swap=False):
        npart = lib.sys.npart
        self.p_ids = np.zeros((npart,))
        list_particle_ids = np.vstack((lib.isys.ityp, lib.isys.iso3)).T[:npart]

        stable = np.isin(lib.isys.ityp[:npart], lib.stables.stabvec)
        self.p_ids = np.array([ ptab.modid2pdg[tuple(l)] 
                            for l in list_particle_ids[stable]])
        self.E = lib.coor.p0[:npart][stable]
        if swap:
            self.pz = -lib.coor.pz[:npart][stable]
        else:
            self.pz = lib.coor.pz[:npart][stable]


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
                print 'stable,', pdg_id
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
