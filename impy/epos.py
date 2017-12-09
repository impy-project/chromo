'''
Created on 03.05.2016

@author: afedynitch
'''

from impy.common import MCRun, MCEvent, EventKinematics, standard_particles
import numpy as np

#=========================================================================
# EPOSMCEvent
#=========================================================================


class EPOSMCEvent(MCEvent):
    def __init__(self, lib, event_config):

        evt = lib.hepevt
        nhep = evt.nhep
        sel = None

        if event_config['charged_only']:
            sel = np.where((evt.isthep[:nhep] == 1) & (
                np.abs(lib.charge_vect(evt.idhep[:nhep])) == 1))
        else:
            sel = np.where(evt.isthep[:evt.nhep] == 1)

        if 'charge_info' in event_config and event_config['charge_info']:
            self.charge = lib.charge_vect(evt.idhep[sel])

        self.p_ids = evt.idhep[sel]
        self.pt2 = evt.phep[0, sel][0]**2 + evt.phep[1, sel][0]**2
        self.pz = evt.phep[2, sel][0]
        self.en = evt.phep[3, sel][0]
        # print zip(self.p_ids, self.pz, self.en)
        # Save for later processing steps
        self.sel = sel
        self.lib = lib

        MCEvent.__init__(self, event_config)


#=========================================================================
# EPOSMCRun
#=========================================================================
class EPOSMCRun(MCRun):
    def __init__(self, libref, **kwargs):

        if not kwargs["event_class"]:
            kwargs["event_class"] = EPOSMCEvent

        MCRun.__init__(self, libref, **kwargs)

        self.generator = "EPOS"
        self.version = "LHC"
        self.needs_init = True
        self.stable_list = []
        self.lib.charge_vect = np.vectorize(self.lib.getcharge)

    def get_sigma_inel(self):
        # Calculate inelastic cross section
        k = self.evkin

        return self.lib.xsection()[1]

    def set_event_kinematics(self, event_kinematics):
        k = event_kinematics
        self.epos_tup = (k.ecm, -1., k.p1pdg, k.p2pdg, k.A1, k.Z1, k.A2, k.Z2)
        self.lib.initeposevt(*self.epos_tup)
        print "Set event kinematics with", k.ecm
        self.event_config['event_kinematics'] = event_kinematics

    def init_generator(self, config, datdir='./iamdata/'):

        # datdir = '../../iamdata/'
        ounit = None
        try:
            from MCVD.management import LogManager  # @UnresolvedImport
            self.log_man = LogManager(config, self)
            self.log_man.create_log('{0}_{1}_{2}'.format(
                self.evkin.A1,
                self.evkin.Z1, self.evkin.p1pdg), '{0}_{1}'.format(
                    self.evkin.A2, self.evkin.Z2), self.evkin.ecm)
            ounit = 66
        except ImportError:
            print(self.class_name +
                  "::init_generator(): Running outside of MCVD," +
                  "the log will be printed to STDOUT.")
            ounit = 6

        if self.debug:
            print self.class_name + "::init_generator(): output_unit set to", \
                ounit
            print self.class_name + "::init_generator(): calling init with:", \
                (1., 1e6, datdir, len(datdir), 1, 2212, 2212, 1, 1, 1, 1, 0, ounit)

        try:
            self.lib.init
        except:
            self.lib.aaset(0)
            self.lib.initializeepos(1., 1e6, datdir,
                                    len(datdir), 1, 2212, 2212, 1, 1, 1, 1, 0,
                                    ounit)
            self.lib.init = True

        # Set default stable
        for pid in [2112, 111, 211, -211, 321, -321, 310, 13, -13]:
            self.set_stable(pid)

        if 'stable' in self.event_config:
            for pdgid in self.event_config['stable']:
                self.set_stable(pdgid)

        # Comprise DPMJET input from the event kinematics object
        self.set_event_kinematics(self.evkin)

    def set_stable(self, pdgid):
        if pdgid not in self.stable_list:
            self.lib.setstable(pdgid)
            self.stable_list.append(pdgid)
        print self.class_name + "::set_stable(): defining ", \
            pdgid, "as stable particle"

    def generate_event(self):

        self.lib.aepos(-1)
        self.lib.afinal()
        self.lib.hepmcstore()

        return False


#=========================================================================
# EPOSMCRun
#=========================================================================


class EPOSCascadeRun():
    def __init__(self,
                 lib_str,
                 label,
                 decay_mode,
                 n_events,
                 fill_subset=False,
                 p_debug=True,
                 inimass=14):

        from ParticleDataTool import DpmJetParticleTable
        exec "import " + lib_str + " as eposlib"
        self.lib = eposlib  # @UndefinedVariable
        self.label = label
        self.nEvents = n_events
        self.spectrum_hists = []

        self.dbg = p_debug
        self.fill_subset = fill_subset
        self.inimass = inimass
        self.ptab = DpmJetParticleTable()
        self.projectiles = [
            'p', 'n', 'K+', 'K-', 'pi+', 'pi-', 'K0L', 'p-bar', 'n-bar',
            'Sigma-', 'Sigma--bar', 'Sigma+', 'Sigma+-bar', 'Xi0', 'Xi0-bar',
            'Xi-', 'Xi--bar', 'Lambda0', 'Lambda0-bar'
        ]
        self.proj_allowed = [
            self.ptab.modname2pdg[p] for p in self.projectiles
        ]
        self.init_generator()
        self.set_stable(decay_mode)

    def get_hadron_air_cs(self, E_lab, projectile_id):
        proj_pdg = self.ptab.modname2pdg[projectile_id]
        # proj_pdg = self.lib.idtrafo("sib","pdg", projectile_id)
        self.epos_tup = (-1., E_lab, proj_pdg, 2212, 1, 1, 14, 7)

        self.lib.initeposevt(*self.epos_tup)
        if self.dbg:
            print(
                "EPOSCascadeRun::get_hadron_air_cs():" +
                "calculating cross-section calculation for {0} @ {1:3.1g} GeV"
            ).format(projectile_id, E_lab)
        # Calculate inelastic cross section
        return self.lib.xsection()[1]

    def init_generator(self):
        from random import randint

        #         datdir = '/lustre/fs17/group/that/af/m2m/iamdata/'
        datdir = '../../iamdata/'
        # Initialize for maximum energy
        self.evkin = EventKinematics(
            plab=1e10,
            p1pdg=2212,
            nuc2_prop=(self.inimass, int(self.inimass / 2)))
        k = self.evkin

        # Set seed of random number generator
        rseed = float(randint(1000000, 10000000))
        print self.__class__.__name__ + '::init_generator(): seed=', rseed
        ounit = 6
        self.lib.aaset(0)
        # Set lab frame "= 2"
        self.lib.initializeepos(rseed, k.ecm, datdir,
                                len(datdir), 2, 2212, 2212, k.A1, k.Z1, k.A2,
                                k.Z2, 0, ounit)

        self.epos_tup = (k.ecm, -1., k.p1pdg, k.p2pdg, k.A1, k.Z1, k.A2, k.Z2)

        self.lib.initeposevt(*self.epos_tup)

    def start(self, projectile, E_lab, sqs, Atarget=0):
        templ = '''
        Model     : {0}
        N_events  : {1}
        Projectile: {2}
        E_lab     : {3:5.3e} GeV
        E_cm      : {3:5.3e} GeV
        Target    : {4}
        '''
        print templ.format(self.label, self.nEvents, projectile, E_lab,
                           Atarget)

        p1pdg = self.ptab.modname2pdg[projectile]
        if self.dbg:
            print 'generating {0} {1} events with {2}'.format(
                self.nEvents, projectile, self.label)

        if p1pdg not in self.proj_allowed:
            raise Exception("EPOSCascadeRun(): Projectile " + int(p1pdg) +
                            " not allowed.")
        if Atarget <= 0:
            self.epos_tup = (-1., E_lab, p1pdg, 2212, 1, 1, 14, 7)
        else:

            self.epos_tup = (-1., E_lab, p1pdg, 2212, 1, 1, Atarget, int(
                Atarget / 2) if Atarget > 1 else 1)

        hist_d = {}
        for hist in self.spectrum_hists:
            hist_d[hist.particle_id] = hist
        ngenerated = self.nEvents
        self.lib.initeposevt(*self.epos_tup)

        for i in xrange(self.nEvents):
            self.lib.aepos(-1)
            self.lib.afinal()
            self.lib.hepmcstore()

            if not (i % 10000) and i and self.dbg:
                print i, "events generated."

            event = EPOSCascadeEvent(self.lib)

            unique_pids = np.unique(event.p_ids)
            if not self.fill_subset:
                [hist_d[pid].fill_event(event) for pid in unique_pids]
            else:
                for pid in unique_pids:
                    if pid in hist_d.keys():
                        hist_d[pid].fill_event(event)
        # Correct for selective filling of histograms
        for hist in self.spectrum_hists:
            hist.n_events_filled = ngenerated

    def set_stable(self, decay_mode):
        from ParticleDataTool import SibyllParticleTable
        if self.dbg:
            print(self.__class__.__name__ + "::set_stable():" +
                  " Setting standard particles stable. " + "Decay mode =",
                  decay_mode)

        # fast-mode particles
        if decay_mode == 0:
            for pdg_id in standard_particles:
                if self.dbg:
                    print 'stable,', pdg_id
                self.lib.setstable(pdg_id)
            return

        # keep muons pions, kaons
        stab = SibyllParticleTable()
        for i in range(4, 5 + 1):
            if self.dbg:
                print 'stable,', stab.modid2pdg[i]
            self.lib.setstable(stab.modid2pdg[i])
        for i in range(7, 18 + 1):
            if self.dbg:
                print 'stable,', stab.modid2pdg[i]
            self.lib.setstable(stab.modid2pdg[i])
        # K0 and K0-bar have to remain unstable to form K0S/L
        # set pi0 unstable

        if decay_mode <= 1:
            return

        # Decay mode 2 for generation of decay spectra (all conventional with
        # lifetime >= K0S
        if self.dbg:
            print(self.__class__.__name__ + "::set_stable(): Setting "
                  "conventional Sigma-, Xi0," +
                  "Xi- and Lambda0 stable (decay mode).")
        for i in range(36, 39 + 1):
            if self.dbg:
                print 'stable,', stab.modid2pdg[i]
            self.lib.setstable(stab.modid2pdg[i])

        if decay_mode <= 2:
            return

        # Conventional mesons and baryons
        # keep eta, eta', rho's, omega, phi, K*
        if self.dbg:
            print(self.__class__.__name__ + "::set_stable(): Setting" +
                  " all conventional stable.")
        # pi0
        if self.dbg:
            print 'stable,', 111
        self.lib.setstable(111)

        for i in range(23, 33 + 1):
            if self.dbg:
                print 'stable,', stab.modid2pdg[i]
            self.lib.setstable(stab.modid2pdg[i])

        # keep SIGMA, XI, LAMBDA
        for i in range(34, 49 + 1):
            if self.dbg:
                print 'stable,', stab.modid2pdg[i]
            self.lib.setstable(stab.modid2pdg[i])

        if decay_mode <= 3:
            return

        if decay_mode <= 4:
            for i in range(59, 99 + 1):
                if i > 60 and i < 71 or i == 82 or \
                        i > 89 and i < 94:
                    continue
                if self.dbg:
                    print 'stable,', stab.modid2pdg[i]
                self.lib.setstable(stab.modid2pdg[i])


#=========================================================================
# EPOSMCEvent
#=========================================================================
class EPOSCascadeEvent():
    def __init__(self, lib):

        evt = lib.hepevt

        sel = np.where(evt.isthep[:evt.nhep] == 1)

        self.p_ids = evt.idhep[sel]
        self.pz = evt.phep[2, sel][0]
        self.E = evt.phep[3, sel][0]
