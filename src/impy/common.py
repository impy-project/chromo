"""The :mod:`common`module contains the classes that expose the
interaction model interface to the front-end and/or the user.

Classes derived from :class:`MCEvent` in (for example :mod:`impy.models.sibyll`)
cast the data from the event generators' particle stacks to numpy arrays.
The basic variables are sufficient to compute all derived attributes,
such as the rapidity :func:`MCEvent.y` or the laboratory momentum fraction
:func:`MCEvent.xlab`.
"""
from abc import ABC, abstractmethod
import numpy as np
from impy import impy_config
from impy.util import info
import typing as _tp
from dataclasses import dataclass


@dataclass
class EventData:
    """Data structure to keep filtered data"""

    nevent: int
    id: np.array
    status: np.array
    charge: np.array
    px: np.array
    py: np.array
    pz: np.array
    en: np.array
    m: np.array
    vx: np.array
    vy: np.array
    vz: np.array
    vt: np.array

    def __getitem__(self, arg):
        """Filter event."""
        if not isinstance(arg, _tp.Iterable):
            raise NotImplementedError
        return EventData(
            self.nevent,
            self.id[arg],
            self.status[arg],
            self.charge[arg],
            self.px[arg],
            self.py[arg],
            self.pz[arg],
            self.en[arg],
            self.m[arg],
            self.vx[arg],
            self.vy[arg],
            self.vz[arg],
            self.vt[arg],
        )

    def __len__(self):
        return len(self.id)


class MCEvent(ABC, EventData):
    """The basis of interaction between user and all the event generators.

    The derived classes are expected to interact with the particle stack
    and derive the base variables from which everything else can be defined.

    Args:
        lib (object)       : Reference to the FORTRAN library in use
        event_kinematics   : Reference to current event kinematics object
        event_frame (str)  : The frame in which the generator returned the variables
    """

    # name of hepevt common block (override in subclass if necessary)
    _hepevt = "hepevt"

    # names of fields of hepevt common block
    # (override in subclass if necessary)
    _nevhep = "nevhep"
    _nhep = "nhep"
    _idhep = "idhep"
    _isthep = "isthep"
    _phep = "phep"
    _vhep = "vhep"
    _jmohep = "jmohep"
    _jdahep = "jdahep"

    # only called once, so setup is allowed to be slow
    def __init__(self, lib, event_kinematics, event_frame):
        self._lib = lib  # used by _charge_init

        # Store for further filtering/access
        self.kin = event_kinematics
        self.event_frame = event_frame

        evt = getattr(lib, self._hepevt)

        npart = getattr(evt, self._nhep)
        sel = slice(None, npart)

        self._parr = getattr(evt, self._phep)[:, sel]
        self._varr = getattr(evt, self._vhep)[:, sel]

        EventData.__init__(
            self,
            getattr(evt, self._nevhep),
            getattr(evt, self._idhep)[sel],
            getattr(evt, self._isthep)[sel],
            self._charge_init(npart),
            self._parr[0],
            self._parr[1],
            self._parr[2],
            self._parr[3],
            self._parr[4],
            self._varr[0],
            self._varr[1],
            self._varr[2],
            self._varr[3],
        )

        if self._jmohep is None:
            self._parents = None
        else:
            self._parents = getattr(evt, self._jmohep)[:, sel]
        if self._jdahep is None:
            self._children = None
        else:
            self._children = getattr(evt, self._jdahep)[:, sel]

        # Apply boosts into frame required by user
        self.kin.apply_boost(self, event_frame, impy_config["user_frame"])
        self.event_frame = impy_config["user_frame"]

    def final_state(self):
        """After calling this method, the variables will only contain
        "stable" final state particles.
        """
        return self.__getitem__(self.status == 1)

    def final_state_charged(self):
        """After calling this method, the variables will only contain
        "stable" and charged final state particles.
        """
        return self.__getitem__((self.status == 1) & (self.charge != 0))

    @abstractmethod
    def _charge_init(self, npart):
        pass

    @property
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
        return self._parents

    @property
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
        return self._children

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
        return 2.0 * self.pz / self.kin.ecm

    @property
    def theta(self):
        """arctan(pt/pz) in rad"""
        return np.arctan2(self.pt, self.pz)

    @property
    def phi(self):
        """arctan(py/px) in rad"""
        return np.arctan2(self.py, self.px)

    @property
    def elab(self):
        """Kinetic energy"""
        kin = self.kin
        if self.event_frame == "laboratory":
            return self.en
        return kin.gamma_cm * self.en + kin.betagamma_z_cm * self.pz

    @property
    def ekin(self):
        """Kinetic energy"""
        return self.elab - self.m

    @property
    def xlab(self):
        """Energy fraction E/E_beam in lab. frame"""
        return self.elab / self.kin.elab

    @property
    def fw(self):
        """I don't remember what this was for..."""
        return self.en / self.kin.pcm


# =========================================================================
# Settings
# =========================================================================
class Settings(ABC):
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
        if self.__class__.__name__ != other_instance.__class__.__name__:
            return True

        other_attr = other_instance.__dict__

        for attr, value in self.__dict__.items():
            if attr == "lib":
                continue
            elif attr not in other_attr.keys():
                return True
            elif value != other_attr[attr]:
                return True
        return False


# =========================================================================
# MCRun
# =========================================================================
class MCRun(ABC):
    #: Prevent creating multiple classes within same python scope
    _is_initialized = []

    def __init__(self, interaction_model_def, settings_dict=dict(), **kwargs):
        import importlib

        # from impy.util import OutputGrabber

        # Import library from library name
        self.lib = importlib.import_module(
            "impy.models." + interaction_model_def.library_name
        )

        # Save definitions from namedtuple into attributes
        self.library_name = interaction_model_def.library_name
        self._event_class = interaction_model_def.EventClass
        self._name = interaction_model_def.name
        self._version = interaction_model_def.version
        self._output_frame = interaction_model_def.output_frame
        if "label" not in kwargs:
            self._label = self._name + " " + self._version
        else:
            self._label = kwargs["label"]

        # Currently initialized event kinematics
        self._curr_event_kin = None

        # Not yet clear how to handle these
        self.setting_dict = settings_dict

        # FORTRAN LUN that keeps logfile handle
        self.output_lun = None

        # Number of generated events so far
        self.nevents = 0

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
    def init_generator(self, event_kinematics, seed="random"):
        """Initializes event generator.

        The maximal energy and particle masses from the event_kinematics
        object define the maximal range, i.e. the energy requested in subsequent
        `_set_event_kinematics` calls should not exceed the one provided here.

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
    def _set_event_kinematics(self, evtkin):
        """Set new combination of energy, momentum, projectile
        and target combination for next event.

        Either, this method defines some derived variables
        that generate_event() can use to generate new events
        without additional arguments, or, it can also set
        internal variables of the model. In both cases the
        important thing is that generate_event remains argument-free.
        """
        # _set_event_kinematics may call library functions,
        # which work correctly only after calling other
        # (initializing) library functions.
        # Therefore _set_event_kinematics should only be called
        # after these calls in init_generator
        pass

    @property
    def event_kinematics(self):
        return self._curr_event_kin

    @event_kinematics.setter
    def event_kinematics(self, evtkin):
        self._set_event_kinematics(evtkin)

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
    def sigma_inel(self, **kwargs):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        pass

    # TODO: Change to generic function for composite target. Make
    # exception for air
    def sigma_inel_air(self, **kwargs):
        """Hadron-air production cross sections according to current
        event setup (energy, projectile).

        Args:
           precision (int): Anything else then 'default' (str) will set
                            the number of MC trails to that number.
        """
        from copy import copy
        from impy.kinematics import EventKinematics

        # Make a backup of the current kinematics config
        prev_kin = copy(self._curr_event_kin)
        frac_air = impy_config["frac_air"]

        cs = 0.0
        for f, iat in frac_air:
            if prev_kin.p1_is_nucleus:
                k = EventKinematics(
                    ecm=prev_kin.ecm,
                    nuc1_prop=(prev_kin.A1, prev_kin.Z2),
                    nuc2_prop=(iat, int(iat / 2)),
                )
            else:
                k = EventKinematics(
                    ecm=prev_kin.ecm,
                    p1pdg=prev_kin.p1pdg,
                    nuc2_prop=(iat, int(iat / 2)),
                )
            self._set_event_kinematics(k)
            cs += f * self.sigma_inel(**kwargs)

        # Restore settings
        self._set_event_kinematics(prev_kin)

        return cs

    @abstractmethod
    def attach_log(self, fname):
        """Routes the output to a file or the stdout."""
        pass

    def _abort_if_already_initialized(self):
        """The first initialization should not be run more than
        once. This method should be called in the beginning of each
        init_generator() implementation.
        """
        message = """
        Don't run initialization multiple times for the same generator. This
        is a limitation of fortran libraries since all symbols are by default
        in global scope. Multiple instances can be created in mupliple threads
        or "python executables" using Pool in multiprocessing etc."""

        assert self.library_name not in self._is_initialized, message
        self._is_initialized.append(self.library_name)

    def _attach_fortran_logfile(self, fname):
        """Chooses a random LUN between 20 - 100 and returns a FORTRAN
        file handle (LUN number) to an open file."""
        import os
        from os import path
        from random import randint

        if path.isfile(fname) and os.stat(fname).st_size != 0:
            raise Exception("Attempts to overwrite log :" + fname)
        elif self.output_lun is not None:
            raise Exception("Log already attached to LUN", self.output_lun)

        path.abspath(fname)
        # Create a random fortran output unit
        self.output_lun = randint(20, 100)
        self.lib.impy_openlogfile(path.abspath(fname), self.output_lun)
        return self.output_lun

    def close_logfile(self):
        """Constructed for closing C++ and FORTRAN log files"""
        if "pythia8" not in self._label:
            self.close_fortran_logfile()
        else:
            # self.close_cc_logfile()
            pass

    def close_fortran_logfile(self):
        """FORTRAN LUN has to be released when finished to flush buffers."""
        if self.output_lun is None:
            info(2, "Output went not to file.")
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

        info(5, "Setting following particles to be stable:", impy_config["stable_list"])

        for pdgid in impy_config["stable_list"]:
            self.set_stable(pdgid)

        if impy_config["pi0_stable"]:
            self.set_stable(111)

    def __call__(self, nevents):
        """Generator function (in python sence)
        which launches the underlying event generator
        and returns its the result (event) as MCEvent object
        """
        retry_on_rejection = impy_config["retry_on_rejection"]
        # Initialize counters to prevent infinite loops in rejections
        ntrials = 0
        nremaining = nevents
        while nremaining > 0:
            if self.generate_event() == 0:
                self.nevents += 1
                yield self._event_class(
                    self.lib, self._curr_event_kin, self._output_frame
                )
                nremaining -= 1
                ntrials += 1
            elif retry_on_rejection:
                info(10, "Rejection occured. Retrying..")
                ntrials += 1
                continue
            elif ntrials > 2 * nevents:
                raise Exception("Things run bad. Check your input.")
            else:
                info(0, "Rejection occured")
