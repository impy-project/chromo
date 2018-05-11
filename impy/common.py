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

# TODO: need full HEPEVT record, something equivalent to
#    int        n;               //!< Number of entries in the event
#    int        ist[n];     //!< Status code
#    int        id [n];     //!< PDG ID
#    int        jmo[n][2];  //!< Pointer to position of 1st and 2nd (or last!) mother
#    int        jda[n][2];  //!< Pointer to position of 1nd and 2nd (or last!) daughter
#    momentum_t p  [n][5];  //!< Momentum: px, py, pz, e, m
#    momentum_t v  [n][4];  //!< Time-space position: x, y, z, t


class MCEvent(object):
    """The basis of interaction between user and all the event generators.

    The derived classes are expected to interact with the particle stack
    and derive the base variables from which everything else can be defined.

    The class attribute __sliced_params__ lists the names to parameters
    which are views on the fortran arrays, and that needs to be changed
    when applying the filter functions.

    Args:
        lib (object)       : Reference to the FORTRAN library in use
        event_kinematics   : Reference to current event kinematics object
        event_frame (str)  : The frame in which the generator returned the variables
        nevent (int)       : Number of current event
        npart (int)        : Number of particles on the stack & length of px, py...
        p_ids (np.array)   : particle ID according to PDG scheme
        status (np.array)  : HEPEVT status flag 1=final state
        px (np.array)      : x-momentum in GeV/c
        py (np.array)      : y-momentum in GeV/c
        pz (np.array)      : z-momentum in GeV/c
        en (np.array)      : Energy in GeV
        m (np.array)       : Mass in GeV
        vx (np.array)      : x-vertex in fm, mm?
        vy (np.array)      : y-vertex in fm, mm?
        vz (np.array)      : z-vertex in fm, mm?
        vt (np.array)      : temporal order of vertex in ps, fs?
    """
    __metaclass__ = ABCMeta
    __sliced_params__ = [
        'p_ids', 'status', 'px', 'py', 'pz', 'en', 'm', 'vx', 'vy', 'vz', 'vt'
    ]

    def __init__(self, lib, event_kinematics, event_frame, nevent, npart,
                 p_ids, status, px, py, pz, en, m, vx, vy, vz, vt):

        # Store the variables for further filtering/access
        self.lib = lib
        self.kin = event_kinematics
        self.event_frame = event_frame

        self.nevent = nevent
        self.npart = npart
        self.p_ids = p_ids
        self.status = status
        self.px = px
        self.py = py
        self.pz = pz
        self.en = en
        self.m = m
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.vt = vt
        
        # Initialize current selection to all entries up to npart
        self.selection = slice(None, self.npart)
        self._apply_slicing()

        # Shall we do this?
        # if impy_config['event_scope'] == 'all':
        #     pass
        # elif impy_config['event_scope'] == 'stable':
        #     self.filter_final_state()
        # elif impy_config['event_scope'] == 'charged':
        #     self.filter_final_state_charged()
        # else:
        #     raise Exception('Unknown event scope')

        # Apply boosts into frame required by user
        self.kin.apply_boost(self, event_frame, impy_config["user_frame"])

    def _apply_slicing(self):
        """Slices/copies the all varaibles according to filter criteria"""
        for var in self.__sliced_params__:
            setattr(self, var, getattr(self, var)[self.selection])

    @abstractmethod
    def filter_final_state(self):
        """After calling this method, the variables will only contain
        "stable" final state particles.
        """
        pass

    @abstractmethod
    def filter_final_state_charged(self):
        """After calling this method, the variables will only contain
        "stable" and charged final state particles.
        """
        pass

    @abstractproperty
    def charge(self):
        """Electrical charge"""
        pass

    @abstractproperty
    def mothers(self, p_idx):
        """Range of indices pointing to mother particles.

        Note::
            This doesn't work if event is fltered for stable or charged
            particles only.

        Args:
            p_idx (int) : Integer index of particle in the variables
        
        Returns:
            (tuple)     : (Index of first mother, index of last mother)
        """
        pass
    
    @abstractproperty
    def daughters(self, p_idx):
        """Range of indices pointing to daughter particles.

        Note::
            This doesn't work if event is fltered for stable or charged
            particles only.

        Args:
            p_idx (int) : Integer index of particle in the variables
        
        Returns:
            (tuple)     : (Index of first daughter, index of last daughter)
        """
        pass
    
    @abstractproperty
    def mothers(self):
        """Mother particles.
        
        Note::
            This doesn't work if event is fltered for stable or charged
            particles only. 
        """
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
        if self.frame == 'laboratory':
            return self.en / kin.elab
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
    def frame(self):
        """Returns frame of the final state"""
        pass

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
        """The method to generate a new event.
        
        Returns:
            (int) : Rejection flag = 0 if everything is ok.
        """
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
        init_generator() implementation.
        """

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
            self.lib.impy_closelogfile(self.output_lun)
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
            if self.generate_event() == 0:
                yield self._event_class(self.lib, self._curr_event_kin,
                                        self.frame)
                nremaining -= 1
                ntrials += 1
            elif retry_on_rejection:
                info(10, 'Rejection occured. Retrying..')
                ntrials += 1
                continue
            elif ntrials > 2 * nevents:
                raise Exception('Things run bad. Check your input.')
            else:
                info(0, 'Rejection occured')
