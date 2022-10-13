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
from impy.util import info, classproperty
from impy.kinematics import EventKinematics
import dataclasses
import copy
import typing as _tp
from contextlib import contextmanager


# Objects of this type contain all default initialization directions
# for each interaction model and create dictionaries that link the
# these settings to either the 'tag' or the 'crmc_id'
@dataclasses.dataclass
class InteractionModelDef:
    tag: str
    name: str
    version: str
    crmc_id: int
    library_name: str
    EventClass: type
    output_frame: str


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
        # used by _charge_init and generator-specific methods
        self._lib = generator._lib

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
            children,
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
        data = generator._lib.crranma4
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
        data = generator._lib.crranma4

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
    _set_final_state_particles_called = False
    _evt_kin = None
    nevents = 0  # number of generated events so far

    def __init__(self, seed, logfname):
        import importlib
        from random import randint

        # from impy.util import OutputGrabber

        self._abort_if_already_initialized()

        assert hasattr(self, "_name")
        assert hasattr(self, "_version")
        assert hasattr(self, "_library_name")
        assert hasattr(self, "_event_class")
        assert hasattr(self, "_output_frame")
        self._lib = importlib.import_module(f"impy.models.{self._library_name}")

        # FORTRAN LUN that keeps logfile handle
        self._output_lun = None

        self._attach_log(logfname)

        if seed is None or seed == "random":
            self._seed = randint(1000000, 10000000)
        elif isinstance(seed, int):
            self._seed = seed
        else:
            raise ValueError(f"Invalid seed {seed}")

        info(3, "Using seed:", self._seed)

    def __call__(self, nevents):
        """Generator function (in python sence)
        which launches the underlying event generator
        and returns its the result (event) as MCEvent object
        """
        assert self._set_final_state_particles_called
        retry_on_rejection = impy_config["retry_on_rejection"]
        # Initialize counters to prevent infinite loops in rejections
        ntrials = 0
        nremaining = nevents
        while nremaining > 0:
            self._update_event_kinematics()
            if self._generate_event() == 0:
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

    @property
    def seed(self):
        return self._seed

    @classproperty
    def name(cls):
        """Event generator name"""
        return cls._name

    @classproperty
    def label(cls):
        """Name and version"""
        return f"{cls._name}-{cls._version}"

    @classproperty
    def pyname(cls):
        """Event generator name as it appears in Python code."""
        return cls.__name__

    @classproperty
    def output_frame(cls):
        """Default frame of the output particle stack."""
        return cls._output_frame

    @classproperty
    def version(cls):
        """Event generator version"""
        return cls._version

    @abstractmethod
    def _generate_event(self):
        """The method to generate a new event.

        Returns:
            (int) : Rejection flag = 0 if everything is ok.
        """
        pass

    @abstractmethod
    def _set_event_kinematics(self, evtkin):
        # Set new combination of energy, momentum, projectile
        # and target combination for next event.

        # Either, this method defines some derived variables
        # that _generate_event() can use to generate new events
        # without additional arguments, or, it can also set
        # internal variables of the model. In both cases the
        # important thing is that _generate_event remains argument-free.

        # This call must not update self.event_kinematics, only the
        # generator-specific variables.
        pass

    def _update_event_kinematics(self):
        if self._evt_kin.composite_target:
            evt_kin = self._evt_kin
            if evt_kin.p2_is_nucleus:
                evt_kin.A2, evt_kin.Z2 = evt_kin.composite_target._get_random_AZ()
                self._set_event_kinematics(evt_kin)

    @property
    def random_state(self):
        return RMMARDState()._record_state(self)

    @random_state.setter
    def random_state(self, rng_state):
        rng_state._restore_state(self)

    @property
    def event_kinematics(self):
        return self._evt_kin

    @event_kinematics.setter
    def event_kinematics(self, evt_kin):
        self._evt_kin = evt_kin
        self._set_event_kinematics(evt_kin)

    def set_stable(self, pdgid, stable=True):
        """Prevent decay of unstable particles

        Args:
            pdgid (int)        : PDG ID of the particle
            stable (bool)      : If `False`, particle is allowed to decay
        """
        # put common code here
        if stable:
            info(5, pdgid, "allowed to decay")
        else:
            info(5, "defining", pdgid, "as stable particle")
        self._set_stable(pdgid, stable)

    def set_unstable(self, pdgid):
        """Convenience funtion for `self.set_stable(..., stable=False)`

        Args:
            pdgid(int)         : PDG ID of the particle
        """
        self.set_stable(pdgid, False)

    def sigma_inel(self, evt_kin=None):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        return self._sigma_inel(self.event_kinematics if evt_kin is None else evt_kin)

    @abstractmethod
    def _sigma_inel(self, evtkin):
        pass

    # TODO: Change to generic function for composite target. Make
    # exception for air
    def sigma_inel_air(self):
        """Hadron-air production cross sections according to current
        event setup (energy, projectile).

        Args:
           precision (int): Anything else then 'default' (str) will set
                            the number of MC trails to that number.
        """

        # Make a backup of the current kinematics config
        frac_air = impy_config["frac_air"]

        cs = 0.0
        for f, iat in frac_air:
            k = EventKinematics(
                ecm=self.event_kinematics.ecm,
                particle1=self.event_kinematics.particle1,
                particle2=(iat, int(iat / 2)),
            )
            cs += f * self._sigma_inel(k)

        return cs

    @abstractmethod
    def _attach_log(self, fname):
        """Routes the output to a file or the stdout."""
        pass

    @abstractmethod
    def _set_stable(self, pidid, stable):
        ...

    def _abort_if_already_initialized(self):
        # The first initialization should not be run more than
        # once.
        message = """
        Don't run initialization multiple times for the same generator. This
        is a limitation of fortran libraries since all symbols are by default
        in global scope. Multiple instances can be created in mupliple threads
        or "python executables" using Pool in multiprocessing etc."""

        assert self._library_name not in self._is_initialized, message
        self._is_initialized.append(self._library_name)

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
        self._output_lun = randint(20, 100)
        self._lib.impy_openlogfile(path.abspath(fname), self.output_lun)
        return self._output_lun

    def _close_logfile(self):
        """Constructed for closing C++ and FORTRAN log files"""
        if "pythia8" not in self._label:
            self._close_fortran_logfile()
        else:
            # self.close_cc_logfile()
            pass

    def _close_fortran_logfile(self):
        """FORTRAN LUN has to be released when finished to flush buffers."""
        if self.output_lun is None:
            info(2, "Output went not to file.")
        else:
            self._lib.impy_closelogfile(self.output_lun)
            self._output_lun = None

    def _set_final_state_particles(self):
        """Defines particles as stable for the default 'tau_stable'
        value in the config."""
        # info(5, 'Setting default particles stable with lifetime <',
        #      impy_config['tau_stable'], 's')

        # for pdgid in make_stable_list(impy_config['tau_stable'], pdata):
        #     self.set_stable(pdgid)

        info(5, "Setting following particles to be stable:", impy_config["stable_list"])

        for pdgid in impy_config["stable_list"]:
            self._set_stable(pdgid, True)

        if impy_config["pi0_stable"]:
            self._set_stable(111, True)

        self._set_final_state_particles_called = True

    @contextmanager
    def _temporary_evt_kin(self, evt_kin):
        prev = copy(self.event_kinematics)
        self.event_kinematics = evt_kin
        yield
        self.event_kinematics = prev
