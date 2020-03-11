'''The :mod:`common`module contains the classes that expose the
interaction model interface to the front-end and/or the user.

Classes derived from :class:`MCEvent` in (for example :mod:`impy.models.sibyll`)
cast the data from the event generators' particle stacks to numpy arrays.
The basic variables are sufficient to compute all derived attributes,
such as the rapidity :func:`MCEvent.y` or the laboratory momentum fraction
:func:`MCEvent.xlab`.
'''
import six
import os
from os.path import abspath
from abc import ABCMeta, abstractmethod, abstractproperty

import numpy as np
from particletools.tables import make_stable_list
from impy import impy_config
from impy.util import info


class MCEvent(object, six.with_metaclass(ABCMeta)):
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
    __sliced_params__ = [
        'p_ids', 'status', 'px', 'py', 'pz', 'en', 'm', 'vx', 'vy', 'vz', 'vt',
        'pem_arr', 'vt_arr'
    ]

    def __init__(self, lib, event_kinematics, event_frame, nevent, npart,
                 p_ids, status, px, py, pz, en, m, vx, vy, vz, vt, pem_arr,
                 vt_arr):
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

        # Full arrays of kinematical vectors
        self.pem_arr = pem_arr  # (px, py, pz, E, m)
        self.vt_arr = vt_arr  # (vx, vy, vz, t)

        # Initialize current selection to all entries up to npart
        if impy_config['pre_slice']:
            info(10, 'Pre-slice enabled.')
            self.selection = slice(None, self.npart)
            self._apply_slicing()
        else: 
            info(10, 'Pre-slice disabled.')
            self.selection = slice(None, None)
        
        # The default slice only cuts limits the view to the array to
        # to the current number of entries
        self._is_filtered = False

        # Apply boosts into frame required by user
        self.kin.apply_boost(self, event_frame, impy_config["user_frame"])
        self.event_frame = impy_config["user_frame"]

    def _apply_slicing(self):
        """Slices/copies the all varaibles according to filter criteria"""
        info(30, 'Slicing attributes.')
        for var in self.__sliced_params__:
            # TODO: AF: Not clear if the kinematical arrays should be
            # exposed to the user and sliced at all. If not remove the braching
            # and delete *_arr from sliced params class variable.
            if var[-4:] == '_arr':
                setattr(self, var, getattr(self, var)[:,self.selection])
            else:
                setattr(self, var, getattr(self, var)[self.selection])
        self._is_filtered = True

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

    @abstractmethod
    def parents(self):
        """Range of indices pointing to mother particles.

        The range of children particles is given by two the two
        indices of the 0th axis.

        Note::
            This property has to raise an exception of filtering
            has been applied, otherwise the indices do not point to
            the right positions on the particle stack.

        Returns:
            (array)     : [2, npart] array
            
        Raises:
            (Exception) : if filtering has been applied.
        """
        if self._is_filtered:
            raise Exception(
                'Parent indices do not point to the' +
                'correct particle indices if slicing/filtering is applied.')

    @abstractmethod
    def children(self):
        """Range of indices pointing to daughter particles.
        
        The range of children particles is given by two the two
        indices of the 0th axis.

        Note::
            This property has to raise an exception of filtering
            has been applied, otherwise the indices do not point to
            the right positions on the particle stack.

        Returns:
            (array)     : [2, npart] array

        Raises:
            (Exception) : if filtering has been applied.
        """
        if self._is_filtered:
            raise Exception(
                'Child indices do not point to the' +
                'correct particle indices if slicing/filtering is applied.')

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
        if self.event_frame == 'laboratory':
            return self.en / kin.elab
        return (kin.gamma_cm * self.en + kin.betagamma_z_cm * self.pz) / kin.elab
    
    @property
    def ekin(self):
        """Kinetic energy"""
        return self.en - self.m

    @property
    def fw(self):
        """I don't remember what this was for..."""
        return self.en / self.kin.pcm


#=========================================================================
# Settings
#=========================================================================
class Settings(six.with_metaclass(ABCMeta)):
    """Custom classes derived from this template allow to set certain low
    level variables in the generators before or after initialization, or for
    each event.

    Note::

        This is only relevant for model developers rather than end users.

    """

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

        for attr, value in six.iteritems(self.__dict__):
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
class MCRun(six.with_metaclass(ABCMeta)):
    #: Prevent creating multiple classes within same python scope
    _is_initialized = False

    def __init__(self, interaction_model_def, settings_dict=dict(), **kwargs):
        import importlib
        from impy.util import OutputGrabber

        # Import library from library name
        self.lib = importlib.import_module(interaction_model_def.library_name)

        # Save definitions from namedtuple into attributes
        self._event_class = interaction_model_def.EventClass
        self._name = interaction_model_def.name
        self._version = interaction_model_def.version
        self._output_frame = interaction_model_def.output_frame
        if 'label' not in kwargs:
            self._label = self._name + " " + self._version
        else:
            self._label = kwargs['label']

        # Flag to control if initialization has been already executed
        self._is_initialized = False

        # Not yet clear how to handle these
        self.setting_dict = settings_dict

        # FORTRAN LUN that keeps logfile handle
        self.output_lun = None

    def __enter__(self):
        """TEMP: It would be good to actually use the with construct to
        open and close logfiles on init."""
        # TODO: this is a bug.
        self.attach_log(None)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """This needs to be tested in more complex scenarios..."""
        self.close_logfile()

    @property
    def label(self):
        """Name + version or custom via keyword arg"""
        return self._label

    @property
    def output_frame(self):
        """Default frame of the output particle stack."""
        return self._output_frame

    @property
    def name(self):
        """Event generator name"""
        return self._name

    @property
    def version(self):
        """Event generator version"""
        return self._version

    @abstractmethod
    def init_generator(self, event_kinematics, seed='random'):
        """Initializes event generator.
        
        The maximal energy and particle masses from the event_kinematics
        object define the maximal range, i.e. the energy requested in subsequent
        `set_event_kinematics` calls should not exceed the one provided here.

        Args:
            event_kinematics (object): maximal energy and masses for subsequent runs
            seed (int)               : random seed, at least 8 digit int

        """ 
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
    def set_stable(self, pdgid, stable=True):
        """Prevent decay of unstable particles
        
        Args:
            pdgid (int)        : PDG ID of the particle
            stable (bool)      : If `False`, particle is allowed to decay
        """
        pass
    
    def set_unstable(self, pdgid):
        """Convenience funtion for `self.set_stable(..., stable=False)`
        
        Args:
            pdgid(int)         : PDG ID of the particle
        """
        self.set_stable(pdgid, stable=False)

    @abstractmethod
    def sigma_inel(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        pass

    @abstractmethod
    def attach_log(self, fname):
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
        import os
        from os import path
        from random import randint

        if path.isfile(fname) and os.stat(fname).st_size != 0:
            raise Exception('Attempts to overwrite log :' + fname)
        elif self.output_lun is not None:
            raise Exception('Log already attached to LUN', self.output_lun)

        path.abspath(fname)
        # Create a random fortran output unit
        self.output_lun = randint(20, 100)
        self.lib.impy_openlogfile(path.abspath(fname), self.output_lun)
        return self.output_lun

    def close_logfile(self):
        """Constructed for closing C++ and FORTRAN log files
        """
        if 'pythia8' not in self._label:
            self.close_fortran_logfile()
        else:
            # self.close_cc_logfile()
            pass

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
        # info(5, 'Setting default particles stable with lifetime <',
        #      impy_config['tau_stable'], 's')

        # for pdgid in make_stable_list(impy_config['tau_stable'], pdata):
        #     self.set_stable(pdgid)

        info(5, 'Setting following particles to be stable:',
             impy_config['stable_list'])

        for pdgid in impy_config['stable_list']:
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
                                        self._output_frame)
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
