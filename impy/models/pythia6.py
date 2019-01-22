'''
Created on 19.01.2015

@author: afedynitch
'''

import numpy as np
from impy.common import MCRun, MCEvent


#=========================================================================
# PhoEvent
#=========================================================================
class PYTHIAEvent(MCEvent):
    def __init__(self, lib, event_config):

        evt = lib.pyjets
        nhep = evt.n
        sel = None

        if event_config['charged_only']:
            lib.pyedit(3)
        else:
            lib.pyedit(2)

        if 'charge_info' in event_config and event_config['charge_info']:
            self.charge = [lib.pychge(evt.k[i, 1]) for i in xrange(nhep)]

        self.p_ids = evt.k[:nhep, 1]
        self.pt2 = evt.p[:nhep, 0]**2 + evt.p[:nhep, 1]**2
        self.pz = evt.p[:nhep, 2]
        self.en = evt.p[:nhep, 3]

        # Save for later processing steps
        self.sel = sel
        self.lib = lib

        MCEvent.__init__(self, event_config)


#=========================================================================
# PhoRun
#=========================================================================
class PYTHIAMCRun(MCRun):
    def __init__(self, libref, **kwargs):
        from ParticleDataTool import PYTHIAParticleData
        if not kwargs["event_class"]:
            kwargs["event_class"] = PYTHIAEvent

        MCRun.__init__(self, libref, **kwargs)

        self.generator = "PYTHIA"
        self.version = "6.4.28"
        self.sigmax = 0.0
        self.ptab = PYTHIAParticleData()
        # The threshold defines the log(E) distance between the initialization
        # energy and current energy.

    def get_sigma_inel(self):
        return self.lib.pyint7.sigt[0, 0, 5]

    def set_stable(self, pdgid):
        kc = self.lib.pycomp(pdgid)
        self.lib.pydat3.mdcy[kc - 1, 0] = 0
        print self.class_name + "::set_stable(): defining ", \
            pdgid, "as stable particle"

    def set_event_kinematics(self, event_kinematics):
        k = event_kinematics
        self.event_config['event_kinematics'] = k
        self.p1_type = self.ptab.name(k.p1pdg)
        self.p2_type = self.ptab.name(k.p2pdg)
        self.ecm = k.ecm

    def init_generator(self, config):
        self.set_event_kinematics(self.evkin)

        try:
            from MCVD.management import LogManager
            # initialize log manager
            self.log_man = LogManager(config, self)
            if self.nondef_stable != None:
                self.log_man.create_log(
                    self.p1_type,
                    self.p2_type,
                    self.ecm,
                    suffix=self.nondef_stable)
            else:
                self.log_man.create_log(self.p1_type, self.p2_type, self.ecm)
            self.lib.pydat1.mstu[10] = 66
        except:
            print self.class_name + "::init_generator(): Running outside of MCVD,", \
                "the log will be printed to STDOUT."
            self.lib.pydat1.mstu[10] = 6

        self.lib.pysubs.msel = 2
        self.lib.pyinit('CMS', self.p1_type, self.p2_type, self.ecm)

        if self.def_settings:
            print self.class_name + \
                "::init_generator(): Using default settings:", \
                self.def_settings.__class__.__name__
            self.def_settings.enable()

        # Set pi0 stable
        self.set_stable(111)

        if 'stable' in self.event_config:
            for pdgid in self.event_config['stable']:
                self.set_stable(pdgid)

    def generate_event(self):
        self.lib.pyevnt()
        return 0


class PYTHIAnewMCRun(PYTHIAMCRun):
    def init_generator(self, config):
        self.lib.pytune(383)
        PYTHIAMCRun.init_generator(self, config)

    def generate_event(self):
        self.lib.pyevnw()
        return 0
