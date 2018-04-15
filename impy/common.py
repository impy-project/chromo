'''
Created on 14.04.2014

@author: afedynitch
'''
import numpy as np
from abc import ABCMeta, abstractmethod, abstractproperty

#===============================================================================
# Standard stable particles for for fast air shower cascade calculation
#===============================================================================
# Particles having an anti-partner
standard_particles = [11, 13, 15, 211, 321, 2212, 2112, 3122, 411, 421, 431]
standard_particles += [-pid for pid in standard_particles]

#unflavored particles
standard_particles += [111, 130, 310, 221, 223, 333]


#===============================================================================
# EventKinematics
#===============================================================================
class EventKinematics():
    def __init__(self,
                 ecm=None,
                 plab=None,
                 p1pdg=None,
                 p2pdg=None,
                 beam=None,
                 nuc1_prop=None,
                 nuc2_prop=None,
                 pdata=None):

        if ecm and plab:
            raise Exception(self.__class__.__name__ +
                            '::init(): Define either ecm or plab')

        if p1pdg and nuc1_prop:
            raise Exception(self.__class__.__name__ +
                            '::init(): Define either particle id or ' +
                            'nuclear protperties for side 1.')

        if p2pdg and nuc2_prop:
            raise Exception(self.__class__.__name__ +
                            '::init(): Define either particle id or ' +
                            'nuclear protperties for side 2.')

        masses = {
            2212: 0.93827,
            2112: 0.93957,
            211: 0.13957,
            321: 0.49368,
            3122: 1.11568,
            22: 0.
        }
        # In case the list is not sufficient load ParticleDataTool
        # on demand
        self.pdata = None
        # Store masses of particles on side 1 and 2
        pmass1, pmass2 = None, None
        # Store average nucleon mass
        mnuc = 0.5 * (masses[2212] + masses[2112])

        if p1pdg:
            try:
                pmass1 = masses[p1pdg]
            except:
                self.load_pdata_tool()
                pmass1 = self.pdata.mass(p1pdg)
            self.A1, self.Z1 = 1, 1
            self.p1pdg = p1pdg
        else:
            pmass1 = mnuc
            self.p1pdg = 2212
            self.A1, self.Z1 = nuc1_prop

        if p2pdg:
            try:
                pmass2 = masses[p2pdg]
            except:
                self.load_pdata_tool()
                pmass2 = self.pdata.mass(p2pdg)

            self.p2pdg = p2pdg
            if p2pdg > 0:
                self.A2, self.Z2 = 1, 1
            else:
                self.A2, self.Z2 = 1, -1
        else:
            pmass2 = mnuc
            self.p2pdg = 2212
            self.A2, self.Z2 = nuc2_prop

        if ecm and not (plab or beam):
            self.ecm = ecm
            self.elab = 0.5 * (ecm**2 - pmass1**2 + pmass2**2) / pmass2
            self.plab = self.e2p(self.elab, pmass1)

        elif plab and not (ecm or beam):
            self.plab = plab
            self.elab = np.sqrt(plab**2 + pmass1**2)
            self.ecm = np.sqrt((self.elab + pmass2)**2 - self.plab**2)

        elif beam and not (ecm or plab):
            self.beam = beam
            p1 = np.array([0, 0, self.e2p(beam[0], pmass1), beam[0]])
            p2 = np.array([0, 0, -self.e2p(beam[1], pmass1), beam[1]])
            s = p1 + p2
            self.ecm = np.sqrt(s[3]**2 - np.sum(s[:3]**2))
            self.elab = 0.5 * (self.ecm**2 - pmass1**2 + pmass2**2) / pmass2
            self.plab = self.e2p(self.elab, pmass1)
            print('EventKinematics: beam, ecms, elab, plab', self.beam,
                  self.ecm, self.elab, self.plab)
        else:
            raise Exception(self.__class__.__name__ +
                            '::init(): Define at least ecm or plab')

        self.pmass1 = pmass1
        self.pmass2 = pmass2

        # compute pcm
        s = self.ecm**2
        self.s = s
        self.pcm = np.sqrt(s**2 - 2 * (s * (pmass1 + pmass2) + pmass1 * pmass2)
                           + pmass1**2 + pmass2**2) / (2 * self.ecm)
        self.gamma_cm = (self.elab + pmass2) / self.ecm
        self.betagamma_cm = self.plab / self.ecm

        self.e_range = []

    def e2p(self, E, m):
        return np.sqrt((E + m) * (E - m))

    def load_pdata_tool(self):
        import ParticleDataTool
        if not self.pdata:
            self.pdata = ParticleDataTool.PYTHIAParticleData()

    def __ne__(self, other):
        return self.__dict__ != other.__dict__

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __hash__(self):
        return hash('_'.join([
            '{0}_{1}'.format(key, self.__dict__[key])
            for key in sorted(self.__dict__.keys())
        ]))

    def __le__(self, other):
        return self.ecm < other.ecm

    def __ge__(self, other):
        return self.ecm > other.ecm

    def __repr__(self):
        ostr = 'Event kinematics:\n'
        ostr += '\tecm      : {0:10.5f}\n'.format(self.ecm)
        ostr += '\tpcm      : {0:10.5f}\n'.format(self.pcm)
        ostr += '\telab     : {0:10.5f}\n'.format(self.elab)
        ostr += '\tplab     : {0:10.5f}\n'.format(self.plab)
        ostr += '\tgamma_cm : {0:10.5f}\n'.format(self.gamma_cm)
        ostr += '\tbgamm_cm : {0:10.5f}\n'.format(self.betagamma_cm)
        ostr += '\tpdgid 1  : {0:10}\n'.format(self.p1pdg)
        ostr += '\tnucprop 1: {0}/{1}\n'.format(self.A1, self.Z1)
        ostr += '\tpdgid 2  : {0:10}\n'.format(self.p2pdg)
        ostr += '\tnucprop 2: {0}/{1}\n'.format(self.A2, self.Z2)

        return ostr


#===============================================================================
# MCEvent
#===============================================================================
class MCEvent(object):
    """The basis of interaction between user and all the event generators.

    The derived classes are expected to interact with the particle stack
    and derive the base variables from which everything else can be defined.

    Args:
        event_config (dict): Parameters passed from above
        lib (object)       : Reference to the FORTRAN library in use
        px (np.array)      : x-momentum in GeV/c
        px (np.array)      : y-momentum in GeV/c
        px (np.array)      : z-momentum in GeV/c
        en (np.array)      : Energy in GeV
        p_ids (np.array)   : particle ID according to PDG scheme
    """
    __metaclass__ = ABCMeta

    def __init__(self, event_config, lib, px, py, pz, en, p_ids, npart):
        self.kin = event_config['event_kinematics']
        self.lib = lib

        self.px = px
        self.py = py
        self.pz = pz
        self.en = en
        self.p_ids = p_ids
        self.npart = npart

    # AF: Since these variables are easy to extract and there is 
    # significant overhead using properties for trivial attributes
    # the definition of these variables has to be enforced by style
    # and abc should be avoided here. However, all derived attributes that
    # involve computations have to be defined either via enforcing
    # abstractproperty or, generic attributes should be in the base class.


    # @abstractproperty
    # def p_ids(self):
    #     """Particle IDs in PDG numbering scheme"""
    #     pass

    # @abstractproperty
    # def px(self):
    #     """x-momentum in GeV/c"""
    #     pass

    # @abstractproperty
    # def py(self):
    #     """y-momentum in GeV/c"""
    #     pass

    # @abstractproperty
    # def pz(self):
    #     """z-momentum in GeV/c"""
    #     pass

    # @abstractproperty
    # def en(self):
    #     """Energy in GeV"""
    #     pass

    @abstractproperty
    def charge(self):
        """Electrical charge"""
        pass


    @property
    def pt(self):
        """Transverse momentum in GeV/c"""
        return np.sqrt(self.px**2 + self.py**2)

    @property
    def pt2(self):
        """Transverse momentum squared in (GeV/c)**2"""
        return self.px**2 + self.py**2

    @property
    def p_tot(self):
        """Total momentum in GeV/c"""
        return np.sqrt(self.pt2 + self.pz**2)

    @property
    def eta(self):
        """Pseudo-rapidity"""
        return np.log((self.p_tot + self.pz) / self.pt)

    @property
    def y(self):
        """True rapidity"""
        return 0.5 * np.log((self.en + self.pz) / (self.en - self.pz))

    @property
    def xf(self):
        """Feynman x_F"""
        return 2. * self.pz / self.kin.ecm

    @property
    def xlab(self):
        """Energy fraction E/E_beam in lab. frame"""
        kin = self.kin
        return (kin.gamma_cm * self.en + kin.betagamma_cm * self.pz) / kin.elab

    @property
    def fw(self):
        """I don't remember what this was for..."""
        return self.en / self.kin.pcm


#=========================================================================
# Settings
#=========================================================================
class Settings():
    __metaclass__ = ABCMeta
    def __init__(self, lib):
        self.lib = lib
        self.override_projectile = None

    def get_label(self):
        return self.__class__.__name__

    @abstractmethod
    def enable(self):
        pass

    @abstractmethod
    def reset(self):
        pass

    def __eq__(self, other_instance):
        return not self.__ne__(other_instance)

    def __ne__(self, other_instance):
        if self.__class__.__name__ != \
                other_instance.__class__.__name__:
            return True

        other_attr = other_instance.__dict__

        for attr, value in self.__dict__.iteritems():
            if attr == 'lib':
                continue
            elif attr not in other_attr.keys():
                return True
            elif value != other_attr[attr]:
                return True
        return False


#=========================================================================
# ScanSettings
#=========================================================================
class ScanSettings(Settings):
    @abstractmethod
    def set_current_value(self, value):
        pass


#=========================================================================
# MCRun
#=========================================================================
class MCRun():
    __metaclass__ = ABCMeta
    def __init__(self, libref, label, n_events, DEBUG, event_class,
                 event_config, event_kinematics, default_settings, settings):
        self.lib = libref
        self.label = label
        self.n_events = n_events
        self.triggers = {}
        self.generator = None
        self.version = None
        self.debug = DEBUG
        self.output_unit = -1
        self.event_class = event_class
        self.event_config = event_config
        self.evkin = event_kinematics
        self.class_name = self.__class__.__name__
        self._is_initialized = False

        def_settings_class, def_settings_args = default_settings
        if def_settings_class:
            self.def_settings = def_settings_class(self.lib, def_settings_args)
        else:
            self.def_settings = Settings(self.lib)

        settings_class, settings_args = settings
        if settings_class:
            self.settings = settings_class(self.lib, settings_args)
            print self.class_name + "::__init__(): Attaching custom settings", \
                  self.settings.get_label()
        else:
            self.settings = Settings(self.lib)

        # Make sure that only charged particles are kept
        if not self.event_config:
            self.event_config = {'charged_only': True}
        if 'charged_only' not in self.event_config:
            self.event_config['charged_only'] = True

        self.nondef_stable = None
        if 'stable' in self.event_config.keys():
            self.nondef_stable = np.sum(
                np.abs(np.array(self.event_config['stable'])))
        self.event_config['event_kinematics'] = self.evkin

    @abstractmethod
    def init_generator(self, config):
        pass

    @abstractmethod
    def generate_event(self):
        pass

    @abstractmethod
    def set_event_kinematics(self, evtkin):
        pass

    @abstractmethod
    def set_stable(self, pdgid):
        pass

    @abstractmethod
    def get_sigma_inel(self):
        pass

    def abort_if_already_initialized(self):
        assert not self._is_initialized
        self._is_initialized = True

    def get_model_label(self):
        return self.label

    def attach_log(self, log_fname):
        pass

    def deattach_log(self):
        pass

    def trigger_event(self, event):
        return [
            trigger for trigger in self.triggers.itervalues()
            if trigger.trigger_event(event)
        ]

    def _init_progress_bar(self, maxval=None):
        if maxval == None:
            maxval = self.n_events + 1
        try:
            from progressbar import ProgressBar, Percentage, Bar, ETA
        except ImportError:
            print "Failed to import 'progressbar' -- disabling progress indicator."
            print "Install the module with 'easy_install progressbar', or",
            print "get it from http://qubit.ic.unicamp.br/~nilton"
            raise Exception("It's easy do do this...")
        self.progressBar = ProgressBar(
            maxval=maxval, widgets=[Percentage(), ' ', Bar(), ' ', ETA()])

    def start(self):
        print "running start()"
        # Look for a definition of energy sweep the event kinematics object
        if self.evkin.e_range:
            self.start_variable()
        else:
            # Start fixed energy run
            self.start_fixed()

    def start_fixed(self):
        print "running start_fixed()"
        if not len(self.triggers):
            raise Exception('No trigger defined in ' + self.class_name + "!")

        print self.class_name + \
        "::start_fixed(): starting generation of", self.n_events, "events"
        print self.evkin
        self.settings.enable()

        self._init_progress_bar()
        self.progressBar.start()

        rej_counter = 0
        for i in xrange(self.n_events):  # @UnusedVariable
            reject = self.generate_event()
            self.progressBar.update(self.progressBar.currval + 1)
            if reject:
                rej_counter += 1
                continue

            event = self.event_class(self.lib, self.event_config)
            fired_triggers = self.trigger_event(event)
            if not fired_triggers:
                continue

            active_histograms = [
                hist
                for trigger in fired_triggers
                for hist in trigger.histogram_list.values()
            ]
            [hist.fill_event(event) for hist in active_histograms]

        # Finalize run
        sigmax = self.get_sigma_inel()
        sigma_gen = sigmax * (float(self.n_events - rej_counter)) / \
            float(self.n_events)

        for trigger in self.triggers.itervalues():
            for hist in trigger.histogram_list.itervalues():
                hist.finalize_run(sigma_gen)

        self.progressBar.finish()
        self.settings.reset()
        try:
            self.lib.pho_event(-2, self.p1, self.p2)[1]
        except:
            pass
        print(r"...completed {0:3.2f}\% of the events have been rejected."
              ).format(100. * float(rej_counter) / float(self.n_events))

        try:
            self.log_man.close_log()
        except AttributeError:
            print 'No logging facility defined.'

    def single_event(self):
        #        print "single_event()"
        if not len(self.triggers):
            raise Exception('No trigger defined in ' + self.class_name + "!")

#        print self.evkin

#        self.settings.enable()

        self.generate_event()
        event = self.event_class(self.lib, self.event_config)

        fired_triggers = self.trigger_event(event)
        if not fired_triggers:
            return

        active_histograms = [
            hist
            for trigger in fired_triggers
            for hist in trigger.histogram_list.values()
        ]
        [hist.fill_event(event) for hist in active_histograms]

        # Finalize run

    #        sigmax = self.get_sigma_inel()
    #        print sigmax

    def start_variable(self):
        if not len(self.triggers):
            raise Exception('No trigger defined in ' + self.class_name + "!")

        print self.class_name + \
        "::start_variable(): starting generation of", self.n_events, "events"
        print self.evkin

        self.settings.enable()

        maxval = self.n_events * len(self.evkin.e_range)
        self._init_progress_bar(maxval)
        self.progressBar.start()
        for evkin in self.evkin.e_range:
            self.set_event_kinematics(evkin)
            rej_counter = 0
            for i in xrange(self.n_events):  # @UnusedVariable
                reject = self.generate_event()
                self.progressBar.update(self.progressBar.currval + 1)
                if reject:
                    rej_counter += 1
                    continue

                event = self.event_class(self.lib, self.event_config)
                fired_triggers = self.trigger_event(event)

                if not fired_triggers:
                    continue

                active_histograms = [
                    hist
                    for trigger in fired_triggers
                    for hist in trigger.histogram_list.values()
                ]
                [
                    hist.fill_event(event, evkin.ecm)
                    for hist in active_histograms
                ]

            # Finalize run
            sigmax = self.get_sigma_inel()
            sigma_gen = sigmax * (float(self.n_events - rej_counter)) / \
                float(self.n_events)

            for trigger in self.triggers.itervalues():
                for hist in trigger.histogram_list.itervalues():
                    hist.finalize_run(sigma_gen, evkin.ecm)
            print "...completed. " + str(100 * float(rej_counter) / float(self.n_events)) + \
                  "% of the events have been rejected."

        self.progressBar.finish()
        self.settings.reset()
        self.log_man.close_log()
        del self.lib
