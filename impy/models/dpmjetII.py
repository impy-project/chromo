"""
Created on 17.03.2014

@author: afedynitch
"""

import numpy as np
from common import MCRun, MCEvent, Settings


# ===============================================================================
# DpmjetIIIMCEvent
# ===============================================================================
class DpmjetIIMCEvent(MCEvent):
    def __init__(self, lib, event_config):

        evt = lib.hkkevt
        npart = evt.nhkk
        sel = None
        if event_config["charged_only"]:
            sel = (evt.isthkk[:npart] == 1) & (lib.extevt.idch[:npart] != 0)
        else:
            sel = evt.isthkk[:npart] == 1

        if "charge_info" in event_config and event_config["charge_info"]:
            self.charge = lib.extevt.idch[sel]

        self.p_ids = evt.idhkk[sel]
        self.pt2 = evt.phkk[0, sel] ** 2 + evt.phkk[1, sel] ** 2
        self.pz = evt.phkk[2, sel]
        self.en = evt.phkk[3, sel]

        MCEvent.__init__(self, event_config)


# ===============================================================================
# DpmjetIIMCRun
# ===============================================================================
class DpmjetIIMCRun(MCRun):
    def __init__(self, libref, **kwargs):

        if not kwargs["event_class"]:
            kwargs["event_class"] = DpmjetIIMCEvent

        MCRun.__init__(self, libref, **kwargs)

        self.generator = "DPMJET"
        self.version = "II.55"
        self.needs_init = True

    def get_sigma_inel(self):
        k = self.evkin
        return self.lib.siinel(self.lib.mcihad(k.p1pdg), 1, k.ecm)

    def set_event_kinematics(self, event_kinematics):
        k = event_kinematics

        self.dpmevt_tup = k.elab, k.A1, k.Z1, k.A2, k.Z2, self.lib.mcihad(k.p1pdg)

        self.event_config["event_kinematics"] = k

    def init_generator(self, config):
        # Comprise DPMJET input from the event kinematics object
        self.set_event_kinematics(self.evkin)
        k = self.evkin

        try:
            from MCVD.management import LogManager  # @UnresolvedImport

            self.log_man = LogManager(config, self)
            self.log_man.create_log(
                "{0}_{1}_{2}".format(k.A1, k.Z1, self.lib.mcihad(k.p1pdg)),
                "{0}_{1}".format(k.A2, k.Z2),
                k.ecm,
            )
        except ImportError:
            print(
                self.class_name + "::init_generator(): Running outside of MCVD,",
                "the log will be printed to STDOUT.",
            )
        except AttributeError:
            print(self.class_name + "::init_generator(): Logging not supported.")

        datdir = "/lustre/fs17/group/that/af/m2m/iamdata/"
        from random import randint

        try:
            if self.init_done:
                raise Exception(
                    self.class_name + "::init_generator(): Init already performed"
                )
        except:
            pass
        seed = randint(1000000, 10000000)
        print(self.__class__.__name__ + "::init_generator(): seed=", seed)
        self.lib.dpmjin(seed, datdir)

        if self.def_settings:
            print(
                self.class_name + "::init_generator(): Using default settings:",
                self.def_settings.__class__.__name__,
            )
            self.def_settings.enable()
        self.init_done = True

        # Set pi0 stable
        self.set_stable(111)

        if "stable" in self.event_config:
            for pdgid in self.event_config["stable"]:
                self.set_stable(pdgid)

    def set_stable(self, pdgid):
        kc = self.lib.pycomp(pdgid)
        self.lib.pydat3.mdcy[kc - 1, 0] = 0
        print(
            self.class_name + "::set_stable(): defining ", pdgid, "as stable particle"
        )

    def generate_event(self):
        self.lib.dpmjet2_event(*self.dpmevt_tup, ievframe=2)
        return 0


# ===============================================================================
# DpmjetIICascadeRun
# ===============================================================================
class DpmjetIICascadeRun:
    def __init__(self, lib_str, label, decay_mode, n_events, datdir=None):
        from ParticleDataTool import DpmJetParticleTable

        exec("import " + lib_str + " as dpmlib")
        self.lib = dpmlib  # @UndefinedVariable
        self.label = label
        self.nEvents = n_events
        self.spectrum_hists = []
        self.projectiles = [
            "p",
            "n",
            "K+",
            "K-",
            "pi+",
            "pi-",
            "K0L",
            "p-bar",
            "n-bar",
            "Sigma-",
            "Sigma--bar",
            "Xi0",
            "Xi0-bar",
            "Xi-",
            "Xi--bar",
            "Lambda0",
            "Lambda0-bar",
        ]

        self.ptab = DpmJetParticleTable()
        if datdir == None:
            datdir = "/lustre/fs17/group/that/af/m2m/iamdata/"
        self.init_generator(datdir)
        self.set_stable(decay_mode)

    def init_generator(self, datdir):
        from random import randint

        self.lib.dpmjin(randint(1000000, 10000000), datdir)

    def set_stable(self, decay_mode=1):
        from ParticleDataTool import SibyllParticleTable

        sibtab = SibyllParticleTable()
        print("DpmjetIICascadeRun::set_stable(): Setting standard particles stable.")
        # Define PYTHIA behavior to follow the switches below
        self.lib.pypars.mstp[40] = 2
        # keep muons pions, kaons
        for i in range(4, 5 + 1):
            kc = self.lib.pycomp(sibtab.modid2pdg[i])
            self.lib.pydat3.mdcy[kc - 1, 0] = 0
        for i in range(7, 18 + 1):
            kc = self.lib.pycomp(sibtab.modid2pdg[i])
            self.lib.pydat3.mdcy[kc - 1, 0] = 0
        # K0 and K0-bar have to remain unstable to form K0S/L

        if decay_mode <= 1:
            return

        # Decay mode 2 for generation of decay spectra (all conventional with
        # lifetime >= K0S
        print(
            "SibRun::set_stable(): Setting conventional Sigma-, Xi0,",
            "Xi- and Lambda0 stable (decay mode).",
        )
        for i in range(36, 39 + 1):
            kc = self.lib.pycomp(sibtab.modid2pdg[i])
            self.lib.pydat3.mdcy[kc - 1, 0] = 0

        if decay_mode <= 2:
            return

        # Conventional mesons and baryons
        # keep eta, eta', rho's, omega, phi, K*
        print("SibRun::set_stable(): Setting all conventional stable.")
        # pi0
        kc = self.lib.pycomp(111)
        self.lib.pydat3.mdcy[kc - 1, 0] = 0
        for i in range(23, 33 + 1):
            kc = self.lib.pycomp(sibtab.modid2pdg[i])
            self.lib.pydat3.mdcy[kc - 1, 0] = 0

        # keep SIGMA, XI, LAMBDA
        for i in range(34, 49 + 1):
            kc = self.lib.pycomp(sibtab.modid2pdg[i])
            self.lib.pydat3.mdcy[kc - 1, 0] = 0

        if decay_mode <= 3:
            return

        # Charmed particles (only for version >= 2.2)
        # keep all charmed
        print("SibRun::set_stable(): Setting all conventional and charmed stable.")
        for i in range(59, 99 + 1):
            try:
                kc = self.lib.pycomp(sibtab.modid2pdg[i])
                self.lib.pydat3.mdcy[kc - 1, 0] = 0
            except:
                pass

    def get_hadron_air_cs(self, E_lab, projectile):
        if projectile == "K+" or projectile == "K-" or projectile == "K0L":
            projectile = 3
        elif projectile == "pi+" or projectile == "pi-":
            projectile = 2
        else:
            projectile = 1
        return self.lib.dpjsig(E_lab, projectile)

    def start(self, projectile, E_lab, sqs, Atarget=0):

        print(self.label, self.nEvents, projectile, "events", E_lab, sqs)
        projectile = self.lib.mcihad(self.ptab.modname2pdg[projectile])
        hist_d = {}
        for hist in self.spectrum_hists:
            hist_d[hist.particle_id] = hist
        avail_pid = list(hist_d.keys())
        for i in range(self.nEvents):
            self.lib.dpmjet2_event(E_lab, projectile)

            if not (i % 10000) and i:
                print(i, "events generated.")
            event = DpmjetIICascadeEvent(self.lib)
            #             try:
            [
                hist_d[pid].fill_event(event)
                for pid in np.unique(event.p_ids)
                if pid in avail_pid
            ]
        #             except KeyError:
        #                 print "DpmjetIICascadeRun::start(): Unknown particle id in result."

        # Correct for selective filling of histograms
        for hist in self.spectrum_hists:
            hist.n_events_filled = self.nEvents


# ===============================================================================
# DpmjetIICascadeEvent
# ===============================================================================
class DpmjetIICascadeEvent:
    def __init__(self, lib):

        npart = lib.hkkevt.nhkk
        sel = lib.hkkevt.isthkk[:npart] == 1

        self.p_ids = lib.hkkevt.idhkk[sel]
        self.E = lib.hkkevt.phkk[3, sel]


# ===============================================================================
# StableCharm
# ===============================================================================
class StableCharm(Settings):
    def __init__(self, lib, mode):
        self.mode = mode
        self.prev_settings = {}
        Settings.__init__(self, lib)

    def enable(self):
        self.lib.dechkk(-1)
        if bool(self.prev_settings):
            raise Exception(
                self.__class__.__name__
                + "::enable(): Settings have been already applied"
            )
        print(
            self.__class__.__name__ + "::enable(): " + "Setting normal D-Mesons stable"
        )
        for pdgid in [421, 411, -411, -421]:
            try:
                kc = self.lib.pycomp(pdgid)
                self.prev_settings[kc - 1] = self.lib.pydat3.mdcy[kc - 1, 0]
                self.lib.pydat3.mdcy[kc - 1, 0] = 0
            except:
                pass

        if self.mode <= 1:
            return
        print(
            self.__class__.__name__
            + "::enable(): "
            + "Setting Ds and D(s) resonances stable"
        )
        for pdgid in [431, -431, 433, -433, 423, 413, -413, -423]:
            try:
                kc = self.lib.pycomp(pdgid)
                self.prev_settings[kc - 1] = self.lib.pydat3.mdcy[kc - 1, 0]
                self.lib.pydat3.mdcy[kc - 1, 0] = 0
            except:
                pass

        if self.mode <= 2:
            return
        print(
            self.__class__.__name__ + "::enable(): " + "Setting LambdaC and etaC stable"
        )
        for pdgid in [4122, -4122, 441]:
            try:
                kc = self.lib.pycomp(pdgid)
                self.prev_settings[kc - 1] = self.lib.pydat3.mdcy[kc - 1, 0]
                self.lib.pydat3.mdcy[kc - 1, 0] = 0
            except:
                pass

    def reset(self):
        print(
            self.__class__.__name__ + "::reset(): " + "Restoring charm decay settings."
        )
        for key, value in list(self.prev_settings.items()):
            self.lib.pydat3.mdcy[key, 0] = value
