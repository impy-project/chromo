'''The :mod:`common`module contains the classes that expose the
interaction model interface to the front-end and/or the user.

Classes derived from :class:`MCEvent` in (for example :mod:`impy.models.sibyll`)
cast the data from the event generators' particle stacks to numpy arrays.
The basic variables are sufficient to compute all derived attributes,
such as the rapidity :func:`MCEvent.y` or the laboratory momentum fraction
:func:`MCEvent.xlab`.
'''
import os
from abc import ABCMeta, abstractmethod, abstractproperty
import numpy as np
import yaml

from particletools.tables import PYTHIAParticleData, make_stable_list

# Globals
root_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..")
impy_config = yaml.load(open(os.path.join(root_dir, 'impy_config.yaml')))
pdata = PYTHIAParticleData(
    cache_file=open(
        os.path.join(root_dir, impy_config["pdata_cachefile"]), 'wb'))

from impy.util import info


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
        npart (np.array)   : Number of particle entries on the stack
    """
    __metaclass__ = ABCMeta

    def __init__(self, event_kinematics, lib, px, py, pz, en, p_ids, npart):
        self.kin = event_kinematics
        self.lib = lib

        self.px = px
        self.py = py
        self.pz = pz
        self.en = en
        self.p_ids = p_ids
        self.npart = npart

    # @Hans: Since these variables are easy to extract, and, there is
    # significant overhead using properties for trivial attributes,
    # the definition of these variables has to be enforced by style
    # and abc should be avoided here. However, all derived attributes that
    # involve computations have to be defined either via enforcing
    # abstractproperty or, defined as generic methods below.
    # You can delete this and stuff below after reading, if you agree.
    # Feel free to extend the docstring clarifying this.

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

    @abstractproperty
    def mass(self):
        """Particle mass in GeV as provided by model."""
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
    """Custom classes derived from this template allow to set certain low
    level variables in the generators before or after initialization, or for
    each event.

    Note::

        This is only relevant for model developers rather than end users.

    """
    __metaclass__ = ABCMeta

    def __init__(self, lib):
        self.lib = lib
        # No idea what this was..
        # self.override_projectile = None

    @property
    def label(self):
        """String of the class name for logging."""
        return self.__class__.__name__

    @abstractmethod
    def enable(self):
        """Code, acting on the FORTRAN library :attr:`self.lib` that
        activates some sort of setting."""
        pass

    @abstractmethod
    def reset(self):
        """Code, acting on the FORTRAN library :attr:`self.lib` that
        removes the effect of the activation. 'Reset to default'"""
        pass

    @abstractmethod
    def set_current_value(self, value):
        """Define if you inted to vary some parameter in between
        events."""
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
# MCRun
#=========================================================================
class MCRun():
    __metaclass__ = ABCMeta

    def __init__(
            self,
            libref,
            label=None,
            settings_dict=dict(),
    ):

        self.lib = libref
        self._label = label
        self._is_initialized = False

        # Not yet clear how to handle these
        self.setting_dict = settings_dict

        # FORTRAN LUN that keeps logfile handle
        self.output_lun = None

    def __enter__(self):
        """TEMP: It would be good to actually use the with construct to
        open and close logfiles on init."""
        self.attach_log()
        return self

    # ...

    def __exit__(self, exc_type, exc_value, traceback):
        """This needs to be tested in more complex scenarios..."""
        self.close_fortran_logfile()

    @abstractproperty
    def name(self):
        """Event generator name"""
        pass

    @abstractproperty
    def version(self):
        """Event generator version"""
        pass

    @property
    def label(self):
        """Any tag or label describing the setup"""
        return self._label

    @abstractmethod
    def init_generator(self):
        pass

    @abstractmethod
    def generate_event(self):
        """The method to generate a new event"""
        pass

    @abstractmethod
    def set_event_kinematics(self, evtkin):
        """Set new combination of energy, momentum, projectile
        and target combination for next event.

        Either, this method defines some derived variables
        that generate_event() can use to generate new events
        without additional arguments, or, it can also set
        internal variables of the model. In both cases the
        important thing is that generate_event remains argument-free.
        """
        pass

    @abstractmethod
    def set_stable(self, pdgid):
        """Prevent decay of unstable particles"""
        pass

    @abstractmethod
    def sigma_inel(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        pass

    @abstractmethod
    def attach_log(self):
        """Routes the output to a file or the stdout."""
        pass

    def _abort_if_already_initialized(self):
        """The first initialization should not be run more than
        once. This method should be called in the beginning of each
        init_generator() implementation."""

        assert not self._is_initialized
        self._is_initialized = True

    def _attach_fortran_logfile(self, fname):
        """Chooses a random LUN between 20 - 100 and returns a FORTRAN
        file handle (LUN number) to an open file."""
        from os import path
        from random import randint

        if path.isfile(fname):
            raise Exception('Attempts to overwrite log :' + fname)
        elif self.output_lun is not None:
            raise Exception('Log already attached to LUN', self.output_lun)
        
        path.abspath(fname)
        # Create a random fortran output unit
        self.output_lun = randint(20, 100)
        self.lib.impy_openlogfile(path.abspath(fname), self.output_lun)
        return self.output_lun

    def close_fortran_logfile(self):
        """FORTRAN LUN has to be released when finished to flush buffers."""
        if self.output_lun is None:
            info(2, 'Output went not to file.')
        else:
            self.lib.impy_closelogfile( self.output_lun)
            self.output_lun = None

    def _define_default_fs_particles(self):
        """Defines particles as stable for the default 'tau_stable'
        value in the config."""
        info(5, 'Setting default particles stable with lifetime <',
             impy_config['tau_stable'], 'ps')

        for pdgid in make_stable_list(impy_config['tau_stable'], pdata):
            self.set_stable(pdgid)
            
        if impy_config['pi0_stable']:
            self.set_stable(111)

    def event_generator(self, event_kinematics, nevents):
        """This is some kind of equivalent to Hans'
        generator concept.

        I don't see a good reason to define them as module
        level functions. Due to this fortran library and double
        initialization issue something has to keep track of ownership
        and history. And classes seem to just this.
        """
        self.set_event_kinematics(event_kinematics)
        retry_on_rejection = impy_config['retry_on_rejection']
        # Initialize counters to prevent infinite loops in rejections
        ntrials = 0
        nremaining = nevents
        while nremaining > 0:
            if self.generate_event():
                yield self._event_class(self.lib, self._curr_event_kin)
                nremaining -= 1
                ntrials += 1
            elif retry_on_rejection:
                ntrials += 1
                continue
            elif ntrials > 2 * nevents:
                raise Exception('Things run bad. Check your input.')
            else:
                info(0, 'Rejection occured')
