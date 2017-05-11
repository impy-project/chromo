'''
Created on 17.03.2014

@author: afedynitch
'''
from InteractionModels.common import MCEvent, MCRun, standard_particles
import numpy as np

#===============================================================================
# QGSJet01Run
#===============================================================================
class QGSJet01Run():
    def __init__(self, lib_str, label, decay_mode, n_events, datdir=None):
        from ParticleDataTool import QGSJetParticleTable
        exec "import " + lib_str + " as qgslib"
        self.lib = qgslib  # @UndefinedVariable
        self.label = label
        self.nEvents = n_events
        self.spectrum_hists = []
        if datdir == None:
            datdir = '/lustre/fs17/group/that/af/m2m/iamdata/'
        self.init_generator(datdir)

        self.ptab = QGSJetParticleTable()
        self.projectiles = {'p':2, 'n':2, 'pi+':1, 'pi-':1, 'K+':3, 'K-':3}
        
    def init_generator(self, datdir):
        from random import randint
        seed = randint(1000000, 10000000)
        print self.__class__.__name__ + '::init_generator(): seed=', seed
        self.lib.cqgsini(seed, datdir)
    
    def get_hadron_air_cs(self, E_lab, projectile=1):
        return self.lib.sectnu(E_lab / 14., 14, 1)
    
    def start(self, projectile, E_lab, sqs, Atarget=14):
        if Atarget == 0:
            Atarget = 14
        templ = '''
        Model     : {0}
        N_events  : {1}
        Projectile: {2}
        E_lab     : {3:5.3e} GeV
        E_cm      : {3:5.3e} GeV
        Target    : {4}
        '''
        print templ.format(self.label, self.nEvents, projectile, 
            E_lab, Atarget)

        projectile = self.projectiles[projectile]
        self.lib.xxaini(E_lab, projectile, 1, Atarget)
        
        for i in xrange(self.nEvents):
            self.lib.psconf()
            
            if not (i % 10000) and i:
                print i, "events generated."
            event = QGSJet01Event(self.lib)
          
            [hist.fill_event(event) for hist in self.spectrum_hists]

        for hist in self.spectrum_hists:
            hist.n_events_filled = self.nEvents
        
#===============================================================================
# QGSJet01Event
#===============================================================================
class QGSJet01Event():
    def __init__(self, lib):
        self.p_ids = lib.area14.ich[:lib.area12.nsp]
        self.E = lib.area14.esp[0, :lib.area12.nsp]


#===============================================================================
# QGSJetIIMCRun
#===============================================================================
class QGSJetIIMCRun(MCRun):
    def __init__(self, libref, **kwargs):

        if not kwargs["event_class"]:
            kwargs["event_class"] = QGSJetIIMCEvent
        
        MCRun.__init__(self, libref, **kwargs)

        self.generator = "QGSJET"
        self.version = "II"
        self.needs_init = True
        
        from ParticleDataTool import QGSJetParticleTable
        self.ptab = QGSJetParticleTable()
        # Let the event object know about the particle table
        self.event_config['ptab'] = self.ptab
        
    
    def set_event_kinematics(self, event_kinematics):
        k = event_kinematics
        
        self.qgsproj = self.ptab.pdg2modid[k.p1pdg]
        
        self.event_config['event_kinematics'] = event_kinematics
        self.lib.qgini(k.elab, self.qgsproj, k.A1, k.A2)
        
    
    def init_generator(self, config): 
        from random import randint
        datdir = '/lustre/fs17/group/that/af/m2m/iamdata/'
            
        self.lib.cqgsini(randint(1000000, 10000000), datdir)
        
        self.set_event_kinematics(self.evkin)
        k = self.evkin
        
        try:              
            from MCVD.management import LogManager
            # initialize log manager
            self.log_man = LogManager(config, self)
            self.log_man.create_log(self.qgsproj, k.A2, k.ecm)
            # initialize PHOJET
        except ImportError:
            print self.class_name + "::init_generator(): Running outside of MCVD,", \
                "the log will be printed to STDOUT."
        except AttributeError:
            print self.class_name + "::init_generator(): Logging not supported."    
        
        if self.def_settings:
            print self.class_name + "::init_generator(): Using default settings:", \
                self.def_settings.__class__.__name__
            self.def_settings.enable()
        
        self.lib.qgini(k.elab, self.qgsproj, k.A1, k.A2)
    
    
    def get_sigma_inel(self):
        k = self.evkin
        projectile = None
        if abs(k.p1pdg) in [321, 310, 130] :
            projectile = 3
        elif abs(k.p1pdg) in [211, 111]:
            projectile = 1
        else:
            projectile = 2
        return self.lib.qgsect(k.elab, projectile, k.A1, k.A2)
    
    
    def generate_event(self):
        self.lib.qgconf()
        

#===============================================================================
# QGSJetIIEvent
#===============================================================================
class QGSJetIIMCEvent(MCEvent):
    def __init__(self, lib, event_config):
        
        evt = lib.qgarr14
        npart = lib.qgarr12.nsp
        k = event_config['event_kinematics']
        self.p_ids = evt.ich[:npart]
        sel = np.where(self.p_ids != 0)
        
        charge = np.array([event_config['ptab'].charge_tab[pi] 
                           for pi in self.p_ids])
        self.p_ids = np.array([event_config['ptab'].modid2pdg[pi] 
                               for pi in self.p_ids])

        if event_config['charged_only']:
            sel = np.where(charge != 0)
        
        if 'charge_info' in event_config and event_config['charge_info']:
            self.charge = charge[sel]
            
        self.p_ids = self.p_ids[sel]
        self.px = evt.esp[2, sel][0]
        self.py = evt.esp[3, sel][0]
        self.pt2 = self.px ** 2 + self.py ** 2
        self.plab = evt.esp[1, sel][0]
        self.elab = evt.esp[0, sel][0]
        self.en = k.gamma_cm * self.elab - k.betagamma_cm * self.plab
        self.pz = -k.betagamma_cm * self.elab + k.gamma_cm * self.plab
        
        MCEvent.__init__(self, event_config)

#===============================================================================
# QGSJetIIRun
#===============================================================================
class QGSJetIIRun():
    def __init__(self, lib_str, label, decay_mode, n_events, datdir=None):
        exec "import " + lib_str + " as qgslib"
        self.lib_str = lib_str
        from ParticleDataTool import QGSJetParticleTable
        self.lib = qgslib  # @UndefinedVariable
        self.label = label
        self.nEvents = n_events
        self.spectrum_hists = []
        if datdir == None:
            datdir = '/lustre/fs17/group/that/af/m2m/iamdata/'
        self.init_generator(datdir)
        
        
        self.ptab = QGSJetParticleTable()
        # Model can only treat 3 classes: pion, nucleon, Kaon
        self.projectiles = {'p':2, 'n':3, 'pi+':1, 'pi-':1, 'K+':4, 'K-':-4, 'K0L':5 }
        
        
    def init_generator(self, datdir): 
        from random import randint
        self.lib.cqgsini(randint(1000000, 10000000), datdir)
    
    def get_hadron_air_cs(self, E_lab, projectile):
        if projectile == 'K+' or projectile == 'K-' or projectile == 'K0L':
            projectile = 3
        elif projectile == 'pi+' or projectile == 'pi-':
            projectile = 1
        else:
            projectile = 2
        return self.lib.qgsect(E_lab, projectile, 1, 14)
    
    def start(self, projectile, E_lab, E_cm, Atarget=14):
        if Atarget == 0:
            Atarget = 14
        templ = '''
        Model     : {0}
        N_events  : {1}
        Projectile: {2}
        E_lab     : {3:5.3e} GeV
        E_cm      : {3:5.3e} GeV
        Target    : {4}
        '''
        print templ.format(self.label, self.nEvents, projectile, 
            E_lab, Atarget)
        
        projectile = self.projectiles[projectile]
        if self.lib_str == 'qgsII03' and E_lab > 1e10:
            for hist in self.spectrum_hists:
                hist.n_events_filled = 1
            print "QGSJet-II-03: max. energy 1e10 GeV exceeded..."
            return
        self.lib.qgini(E_lab, projectile, 1, Atarget)
        
        for i in xrange(self.nEvents):
            self.lib.qgconf()
            if not (i % 10000):
                if i: 
                    print i, "events generated."
            event = QGSJetIIEvent(self.lib)
            [hist.fill_event(event) for hist in self.spectrum_hists]

        for hist in self.spectrum_hists:
            hist.n_events_filled = self.nEvents

#===============================================================================
# QGSJetIIEvent
#===============================================================================
class QGSJetIIEvent():
    def __init__(self, lib):
        self.p_ids = lib.qgarr14.ich[:lib.qgarr12.nsp]
        self.E = lib.qgarr14.esp[0, :lib.qgarr12.nsp]
