'''
Created on 15.05.2012

@author: afedynitch
'''
import numpy as np
from impy.common import MCRun, MCEvent

class PhoEvent(MCEvent):
    """Wrapper class around (pure) PHOJET particle stack."""
    def __init__(self, lib, event_config):

        # The stack common block
        evt = lib.poevt1
        # Number of entries on stack
        npart = evt.nhep
        # Filter stack for charged particles if selected
        if event_config['charged_only']:
            sel = np.where((evt.isthep[:npart] == 1) & \
                (np.abs(lib.poevt2.icolor[0, :npart]) == 3))
        else:
            sel = np.where(evt.isthep[:npart] == 1)

        # Save selector for implementation of on-demand properties
        self.sel = sel

        MCEvent.__init__(
            self,
            event_config,
            lib,
            px=evt.phep[0, sel][0],
            py=evt.phep[1, sel][0],
            pz=evt.phep[2, sel][0],
            en=evt.phep[3, sel][0],
            p_ids=evt.idhep[sel],
            npart=npart)

    @property
    def charge(self):
        return self.lib.poevt2.icolor[0, self.sel[0]] / 3

    def gen_cut_info(self):
        """Init variables tracking the number of soft and hard cuts"""

        self._mpi = self.lib.podebg.kspom + self.lib.podebg.khpom
        self._kspom = self.lib.podebg.kspom
        self._khpom = self.lib.podebg.khpom
        self._ksoft = self.lib.podebg.ksoft
        self._khard = self.lib.podebg.khard

    @property
    def mpi(self):
        """Total number of cuts"""
        if not hasattr(self, 'mpi'):
            self.gen_cut_info()
            return self._mpi
        return self._mpi

    @property
    def kspom(self):
        """Total number of soft cuts"""
        if not hasattr(self, 'kspom'):
            self.gen_cut_info()
            return self._kspom
        return self._kspom

    @property
    def khpom(self):
        """Total number of hard cuts"""
        if not hasattr(self, 'khpom'):
            self.gen_cut_info()
            return self._khpom
        return self._khpom

    @property
    def ksoft(self):
        """Total number of realized soft cuts"""
        if not hasattr(self, 'ksoft'):
            self.gen_cut_info()
            return self._ksoft
        return self._ksoft

    @property
    def khard(self):
        """Total number of realized hard cuts"""
        if not hasattr(self, 'khard'):
            self.gen_cut_info()
            return self._khard
        return self._khard

    @property
    def mother(self):
        """The mother particle IDs.

        Mother particles are those with either decay or
        hadronize to the final state particle.
        """
        return self.lib.poevt1.jmohep[0, self.sel]

    @property
    def all_p_ids(self):
        """Unfiltered access to all particle IDs.

        Those include the initial state, quarks, gluons, FSR, ISR, etc.
        """
        return self.lib.poevt1.idhep[:self.npart]

    def elastic_t(self):
        """Squared momentum transfer t for elastic interaction.

        This only makes sense if the interaction is indeed elastic
        and only 4 particles are on the stack. The initial 2 and
        the final 2. Handle with care!!
        """
        return np.sum(
            np.square(lib.poevt1.phep[0:4, 0] - lib.poevt1.phep[0:4, 5]))


# No idea, yet, what to do with this one. It uses some link to
# the fastjet library to do more fancy computations

# class PhoEventJets(PhoEvent):
#     def find_jets(self, dR, psrapcut):
#         # charged anti-kt jets
#         lib = self.lib
#         lib.pho_findjets(1, 1, dR, psrapcut)
#         njets = lib.pojets.njets

#         self.njets = njets
#         self.j_E = lib.pojets.jets[3, :njets]
#         self.j_p_t = np.sqrt(lib.pojets.jets[0, :njets]**2 +
#                              lib.pojets.jets[1, :njets]**2)
#         self.j_pz = lib.pojets.jets[1, :njets]
#         self.j_p_tot = np.sqrt(lib.pojets.jets[0, :njets]**2 +
#                                lib.pojets.jets[1, :njets]**2 +
#                                lib.pojets.jets[2, :njets]**2)

#         self.j_eta = np.log((self.j_p_tot + self.j_pz) / self.j_p_t)
#         self.j_y = 0.5 * \
#             np.log((self.j_E + self.j_pz) / (self.j_E - self.j_pz))


class PHOJETRun(MCRun):
    def __init__(self, libref, **kwargs):

        if not kwargs["event_class"]:
            kwargs["event_class"] = PhoEvent

        MCRun.__init__(self, libref, **kwargs)

        self.generator = "PHOJET"
        self.version = "1.12-35"
        self.sigmax = 0.0
        # The threshold defines the log(E) distance between the initialization
        # energy and current energy.

    def get_sigma_inel(self):
        print "PHOJET workaround for cross-section"
        self.lib.pho_event(1, self.p1, self.p2)[1]
        return self.lib.powght.siggen[3]

    def set_stable(self, pdgid):
        if abs(pdgid) == 2212:
            return
        kc = self.lib.pycomp(pdgid)
        self.lib.pydat3.mdcy[kc - 1, 0] = 0
        print self.class_name + "::set_stable(): defining ", \
            pdgid, "as stable particle"

    def set_event_kinematics(self, event_kinematics):
        k = event_kinematics
        self.event_config['event_kinematics'] = k
        self.p1_type, self.p2_type, self.ecm, self.pcm = \
                            k.p1pdg, k.p2pdg, k.ecm, k.pcm
        if self.def_settings.override_projectile != None:
            print 'Overriding projectile', self.p1_type, self.def_settings.override_projectile
            self.lib.pho_setpar(1, self.def_settings.override_projectile, 0,
                                0.0)
        else:
            self.lib.pho_setpar(1, self.p1_type, 0, 0.0)
        self.lib.pho_setpar(2, self.p2_type, 0, 0.0)
        p1, p2 = np.array(
            np.zeros(4), dtype='d'), np.array(
                np.zeros(4), dtype='d')
        p1[0] = 0.0
        p1[1] = 0.0
        p1[2] = k.pcm
        p1[3] = k.ecm / 2.
        p2[0] = 0.0
        p2[1] = 0.0
        p2[2] = -k.pcm
        p2[3] = k.ecm / 2.

        self.p1, self.p2 = p1, p2

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

            rejection = self.lib.pho_init(
                -2 if self.lib.__name__.find('dpmjet') != -1 else -1, 66)
            self.lib.pydat1.mstu[10] = 66
            # initialize PHOJET
        except ImportError:
            print self.class_name + "::init_generator(): Running outside of MCVD,", \
                "the log will be printed to STDOUT."
            rejection = self.lib.pho_init(-1, 6)
            self.lib.pydat1.mstu[10] = 6
        except (TypeError, AttributeError):
            rejection = self.lib.pho_init(-1)

        if rejection:
            raise Exception(self.class_name + "::init_generator():" +
                            "PHOJET rejected the initialization!")

        process_switch = self.lib.poprcs.ipron
        # non-diffractive elastic scattering (1 - on, 0 - off)
        process_switch[0, 0] = 1
        # elastic scattering
        process_switch[1, 0] = 0
        # quasi-elastic scattering (for incoming photons only)
        process_switch[2, 0] = 1
        # central diffration (double-pomeron scattering)
        process_switch[3, 0] = 1
        # particle 1 single diff. dissociation
        process_switch[4, 0] = 1
        # particle 2 single diff. dissociation
        process_switch[5, 0] = 1
        # double diff. dissociation
        process_switch[6, 0] = 1
        # direct photon interaction (for incoming photons only)
        process_switch[7, 0] = 1

        self.set_event_kinematics(self.evkin)

        rejection = self.lib.pho_event(-1, self.p1, self.p2)[1]

        if rejection:
            raise Exception(self.class_name + "::init_generator():" +
                            "PHOJET rejected the event initialization!")

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
        return self.lib.pho_event(1, self.p1, self.p2)[1]


#=========================================================================
# PhoEventViewer
#=========================================================================
class PhoEventViewer(PhoMCRun):
    def init_generator(self):

        rejection = self.lib.pho_init(-1, 6)
        if rejection:
            raise Exception("PHOJET rejected the initialization!")

        process_switch = self.lib.poprcs.ipron
        # non-diffractive elastic scattering (1 - on, 0 - off)
        process_switch[0, 0] = 1
        # elastic scattering
        process_switch[1, 0] = 0
        # quasi-elastic scattering (for incoming photons only)
        process_switch[2, 0] = 1
        # central diffration (double-pomeron scattering)
        process_switch[3, 0] = 1
        # particle 1 single diff. dissociation
        process_switch[4, 0] = 1
        # particle 2 single diff. dissociation
        process_switch[5, 0] = 1
        # double diff. dissociation
        process_switch[6, 0] = 1
        # direct photon interaction (for incoming photons only)
        process_switch[7, 0] = 1

    def init_event(self):
        reject = 0
        self.p1, self.p2 = self.init_beam_setup(self.sqs, self.p1_type,
                                                self.p2_type)
        # initial beam configuration
        self.settings.enable()
        print "Initializing PHOJET with", self.p1_type, self.p2_type, self.sqs
        self.sigmax, self.reject = self.lib.pho_event(-1, self.p1, self.p2,
                                                      self.sigmax, reject)
        if reject:
            raise Exception("Warning: PHOJET rejected the beam configuration.")
        print self.label, ": Generation 1 event for", self.p1_type, self.p2_type, self.sqs

    def next_event(self):
        cursigma = 0.
        reject = 0
        cursigma, reject = self.lib.pho_event(1, self.p1, self.p2, cursigma,
                                              reject)
        if reject != 0:
            print "Warning: PHOJET rejected the event configuration."
            while reject != 0:
                cursigma, reject = self.lib.pho_event(1, self.p1, self.p2,
                                                      cursigma, reject)
        self.event = self.event_class(self.lib)
