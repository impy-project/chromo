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
from impy.kinematics import EventKinematics
import dataclasses
import copy
import typing as _tp


@dataclasses.dataclass
class EventData:
    """
    Data structure to keep filtered data.

    Unlike the original MCEvent, this structure is picklable.

    generator: (str, str)
        Info about the generator, its name and version.
    kin: EventKinematics
        Info about initial state.
    frame: str
        Current event frame (Lorentz frame).
    nevent: int
        Which event in the sequence.
    pid: 1D array of int
        PDG IDs of the particles.
    status: 1D array of int
        Status code of the particles (1 = final).
    charge: 1D array of double
        Charge in units of elementary charge. Fractional for quarks.
    px: 1D array of double
        X coordinate of momentum in GeV/c.
    py: 1D array of double
        Y coordinate of momentum in GeV/c.
    pz: 1D array of double
        Z coordinate of momentum in GeV/c.
    en: 1D array of double
        Energy in GeV.
    m: 1D array of double
        Generated mass in GeV/c^2.
    vx: 1D array of double
        X coordinate of particle in mm.
    vy: 1D array of double
        Y coordinate of particle in mm.
    vz: 1D array of double
        Z coordinate of particle in mm.
    vt: 1D array of double
        Time of particle in s.
    parents: 2D array of int or None
        This array has shape (N, 2) for N particles. The two numbers
        are a range of indices pointing to parent particle(s). This
        array is not available for all generators and it is not
        available for filtered events. In those cases, parents is None.
    children: 2D array of int or None
        Same as parents.
    """

    generator: _tp.Tuple[str, str]
    kin: EventKinematics
    frame: str
    nevent: int
    pid: np.ndarray
    status: np.ndarray
    charge: np.ndarray
    px: np.ndarray
    py: np.ndarray
    pz: np.ndarray
    en: np.ndarray
    m: np.ndarray
    vx: np.ndarray
    vy: np.ndarray
    vz: np.ndarray
    vt: np.ndarray
    parents: _tp.Optional[np.ndarray]
    children: _tp.Optional[np.ndarray]

    def __getitem__(self, arg):
        """
        Return masked event.

        This may return a copy if the result cannot represented as a view.
        """
        return EventData(
            self.generator,
            self.kin,
            self.frame,
            self.nevent,
            self.pid[arg],
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
            None,
            None,
        )

    def __len__(self):
        """Return number of particles."""
        return len(self.pid)

    def __eq__(self, other):
        """
        Return true if events are equal.

        This is mostly useful for debugging and in unit tests.
        """
        return (
            self.kin == other.kin
            and self.frame == other.frame
            and self.nevent == other.nevent
            and np.all(
                (self.pid == other.pid)
                & (self.status == other.status)
                & (self.charge == other.charge)
                & (self.px == other.px)
                & (self.py == other.py)
                & (self.pz == other.pz)
                & (self.en == other.en)
                & (self.m == other.m)
                & (self.vx == other.vx)
                & (self.vy == other.vy)
                & (self.vz == other.vz)
                & (self.vt == other.vt)
                & (self.en == other.en)
            )
            and (self.parents is None or np.all(self.parents == other.parents))
            and (self.children is None or np.all(self.children == other.children))
        )

    def copy(self):
        """
        Return event copy.
        """
        # this should be implemented with the help of copy
        copies = []
        for obj in dataclasses.astuple(self):
            copies.append(obj.copy() if hasattr(obj, "copy") else obj)

        return EventData(*copies)

    def final_state(self):
        """
        Return filtered event with only final state particles.

        The final state is generator-specific, but usually contains only
        long-lived particles.
        """
        return self[self.status == 1]

    def final_state_charged(self):
        """
        Return filtered event with only charged final state particles.

        The final state is generator-specific, see :meth:`final_state`.
        Additionally, this selects only charged particles which can be
        seen by a tracking detector.
        """
        return self[(self.status == 1) & (self.charge != 0)]

    @property
    def pt(self):
        """Return transverse momentum in GeV/c."""
        return np.sqrt(self.pt2)

    @property
    def pt2(self):
        """Return transverse momentum squared in (GeV/c)**2."""
        return self.px**2 + self.py**2

    @property
    def p_tot(self):
        """Return total momentum in GeV/c."""
        return np.sqrt(self.pt2 + self.pz**2)

    @property
    def eta(self):
        """Return pseudorapidity."""
        return np.log((self.p_tot + self.pz) / self.pt)

    @property
    def y(self):
        """Return rapidity."""
        return 0.5 * np.log((self.en + self.pz) / (self.en - self.pz))

    @property
    def xf(self):
        """Return Feynman x_F."""
        return 2.0 * self.pz / self.kin.ecm

    @property
    def theta(self):
        """Return angle to beam axis (z-axis) in rad."""
        return np.arctan2(self.pt, self.pz)

    @property
    def phi(self):
        """Return polar angle around beam axis (z-axis) in rad."""
        return np.arctan2(self.py, self.px)

    @property
    def elab(self):
        """Return kinetic energy in laboratory frame."""
        kin = self.kin
        if self.frame == "laboratory":
            return self.en
        return kin.gamma_cm * self.en + kin.betagamma_z_cm * self.pz

    @property
    def ekin(self):
        """Return kinetic energy in current frame."""
        return self.elab - self.m

    @property
    def xlab(self):
        """Return energy fraction of beam in laboratory frame."""
        return self.elab / self.kin.elab

    # @property
    # def fw(self):
    #     """I don't remember what this was for..."""
    #     return self.en / self.kin.pcm

    def to_hepmc3(self, genevent=None):
        """
        Convert event to HepMC3 GenEvent.

        After conversion, it is possible to traverse and draw the particle history.
        This requires the optional pyhepmc library.

        Parameters
        ----------
        genevent: GenEvent or None, optional
            If a genevent is passed, its content is replaced with this event
            information. If genevent is None (default), a new genevent is created.
        """
        import pyhepmc  # delay import

        if genevent is None:
            genevent = pyhepmc.GenEvent()

        genevent.from_hepevt(
            event_number=self.nevent,
            px=self.px,
            py=self.py,
            pz=self.pz,
            en=self.en,
            m=self.m,
            pid=self.pid,
            status=self.status,
            parents=self.parents,
            children=self.children,
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            vt=self.vt,
        )

        genevent.run_info = pyhepmc.GenRunInfo()
        genevent.run_info.tools = [self.generator + ("",)]
        return genevent

    # if all required packages are available, add extra
    # method to draw event in Jupyter
    try:
        from pyhepmc import GenEvent

        if hasattr(GenEvent, "_repr_html_"):

            def _repr_html_(self):
                return self.to_hepmc3()._repr_html_()

    except ModuleNotFoundError:
        pass


class MCEvent(EventData, ABC):
    """
    The base class for interaction between user and all event generators.

    Derived classes override base methods, if interaction with the particle stack
    of the generator requires special treatment.
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

    def __init__(self, generator):
        """
        Parameters
        ----------
        generator:
            Generator instance.
        """
        self._lib = generator.lib  # used by _charge_init and generator-specific methods

        evt = getattr(self._lib, self._hepevt)

        npart = getattr(evt, self._nhep)
        sel = slice(None, npart)

        phep = getattr(evt, self._phep)[:, sel]
        vhep = getattr(evt, self._vhep)[:, sel]

        parents = getattr(evt, self._jmohep).T[sel] if self._jmohep else None
        children = getattr(evt, self._jdahep).T[sel] if self._jdahep else None

        user_frame = impy_config["user_frame"]  # TODO take from event_kinematics

        EventData.__init__(
            self,
            (generator.name, generator.version),
            generator.event_kinematics,
            user_frame,
            int(getattr(evt, self._nevhep)),
            getattr(evt, self._idhep)[sel],
            getattr(evt, self._isthep)[sel],
            self._charge_init(npart),
            *phep,
            *vhep,
            parents,
            children
        )

        # Apply boosts into frame required by user
        self.kin.apply_boost(self, generator._output_frame, user_frame)

    @abstractmethod
    def _charge_init(self, npart):
        # override this in derived to get charge info
        ...

    # MCEvent is not pickleable, but EventData is. For convenience, we
    # make it so that MCEvent can be saved and is restored as EventData.
    def __getstate__(self):
        # save only EventData sub-state
        return {k: getattr(self, k) for k in self.__dataclass_fields__}

    def __getnewargs__(self):
        # upon unpickling, create EventData object instead of MCEvent object
        return (EventData,)


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


@dataclasses.dataclass
class RMMARDState:
    _c_number: np.ndarray = None
    _u_array: np.ndarray = None
    _u_i: np.ndarray = None
    _u_j: np.ndarray = None
    _seed: np.ndarray = None
    _counter: np.ndarray = None
    _big_counter: np.ndarray = None
    _sequence_number: np.ndarray = None

    def _record_state(self, generator):
        data = generator.lib.crranma4
        self._c_number = data.c
        self._u_array = data.u
        self._u_i = data.i97
        self._u_j = data.j97
        self._seed = data.ijkl
        self._counter = data.ntot
        self._big_counter = data.ntot2
        self._sequence_number = data.jseq
        return self

    def _restore_state(self, generator):
        data = generator.lib.crranma4

        data.c = self._c_number
        data.u = self._u_array
        data.i97 = self._u_i
        data.j97 = self._u_j
        data.ijkl = self._seed
        data.ntot = self._counter
        data.ntot2 = self._big_counter
        data.jseq = self._sequence_number
        return self

    def __eq__(self, other: object) -> bool:
        if other.__class__ is not self.__class__:
            return NotImplemented

        return (
            np.array_equal(self._c_number, other._c_number)
            and np.array_equal(self._u_array, other._u_array)
            and np.array_equal(self._u_i, other._u_i)
            and np.array_equal(self._u_j, other._u_j)
            and np.array_equal(self._seed, other._seed)
            and np.array_equal(self._counter, other._counter)
            and np.array_equal(self._big_counter, other._big_counter)
            and np.array_equal(self._sequence_number, other._sequence_number)
        )

    def copy(self):
        """
        Return generator copy.
        """
        return copy.deepcopy(self)  # this uses setstate, getstate

    @property
    def sequence(self):
        return self._sequence_number

    @property
    def counter(self):
        return self._counter[self.sequence - 1]

    @property
    def seed(self):
        return self._seed[self.sequence - 1]


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

    def _update_event_kinematics(self):
        if self._curr_event_kin.composite_target:
            ekin = self._curr_event_kin
            if ekin.p2_is_nucleus:
                ekin.A2, ekin.Z2 = ekin.composite_target._get_random_AZ()
                self._set_event_kinematics(ekin)

    @property
    def random_state(self):
        return RMMARDState()._record_state(self)

    @random_state.setter
    def random_state(self, rng_state):
        rng_state._restore_state(self)

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
            self._update_event_kinematics()
            if self.generate_event() == 0:
                self.nevents += 1
                yield self._event_class(self)
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
