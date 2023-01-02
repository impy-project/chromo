"""The :mod:`common`module contains the classes that expose the
interaction model interface to the front-end and/or the user.

Classes derived from :class:`MCEvent` in (for example :mod:`impy.models.sibyll`)
cast the data from the event generators' particle stacks to numpy arrays.
The basic variables are sufficient to compute all derived attributes,
such as the rapidity :func:`MCEvent.y` or the laboratory momentum fraction
:func:`MCEvent.xlab`.
"""
from __future__ import annotations
from abc import ABC, abstractmethod
import numpy as np
from impy.util import (
    classproperty,
    select_parents,
    naneq,
    pdg2name,
    process_particle,
    Nuclei,
)
from impy.constants import (
    quarks_and_diquarks_and_gluons,
    long_lived,
    standard_projectiles,
    GeV,
)
from impy.kinematics import EventKinematics, CompositeTarget
import dataclasses
import copy
from typing import Tuple, Optional
import warnings
from particle import Particle


# Do we need EventData.n_spectators in addition to EventData.n_wounded?
# n_spectators can be computed from n_wounded as number_of_nucleons - n_wounded
# If we want this, it should be computed dynamically via a property.
@dataclasses.dataclass
class CrossSectionData:
    """Information of cross-sections returned by the generator.

    The class is driven by the amount of detail that Pythia-8 provides.
    Most other generators do not fill all attributes. When information is missing,
    the attributes contain NaN. The fields are intentionally redundant, since
    some generators may provide only the total cross-section or only the inelastic
    cross-section.

    All cross-sections are in millibarn.

    Attributes
    ----------
    total : float
        Total cross-section (elastic + inelastic).
    inelastic : float
        Inelastic cross-section. Includes diffractive cross-sections.
    elastic : float
        Cross-section for pure elastic scattering of incoming particle
        without any new particle generation.
    diffractive_xb : float
        Single diffractive cross-section. Particle 2 remains intact,
        particles are produced around Particle 1.
    diffractive_ax : float
        Single diffractive cross-section. Particle 1 remains intact,
        particles are produced around Particle 2.
    diffractive_xx : float
        Double diffractive cross-section. Central rapidity gap, but
        paricles are producted around both beam particles.
    diffractive_axb : float
        Central diffractive cross-section (e.g.
        pomeron-pomeron interaction.)
    """

    total: float = np.nan
    inelastic: float = np.nan
    elastic: float = np.nan
    diffractive_xb: float = np.nan
    diffractive_ax: float = np.nan
    diffractive_xx: float = np.nan
    diffractive_axb: float = np.nan

    @property
    def non_diffractive(self):
        return (
            self.inelastic
            - self.diffractive_xb
            - self.diffractive_ax
            - self.diffractive_xx
            - self.diffractive_axb
        )

    def __eq__(self, other):
        at = dataclasses.astuple(self)
        bt = dataclasses.astuple(other)
        return all(naneq(a, b) for (a, b) in zip(at, bt))

    def __ne__(self, other):
        return not self == other


# Do we need EventData.n_spectators in addition to EventData.n_wounded?
# n_spectators can be computed from n_wounded as number_of_nucleons - n_wounded
# If we want this, it should be computed dynamically via a property.


@dataclasses.dataclass
class EventData:
    """
    Data structure to keep filtered data.

    Unlike the original MCEvent, this structure is picklable.

    generator: (str, str)
        Info about the generator, its name and version.
    kin: EventKinematics
        Info about initial state.
    nevent: int
        Which event in the sequence.
    impact_parameter: float
        Impact parameter for nuclear collisions in mm.
    n_wounded: (int, int)
        Number of wounded nucleons on sides A and B.
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

    generator: Tuple[str, str]
    kin: EventKinematics
    nevent: int
    impact_parameter: float
    n_wounded: Tuple[int, int]
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
    parents: Optional[np.ndarray]
    children: Optional[np.ndarray]

    def __getitem__(self, arg):
        """
        Return masked event.

        This may return a copy if the result cannot represented as a view.
        """
        return EventData(
            self.generator,
            self.kin,
            self.nevent,
            self.impact_parameter,
            self.n_wounded,
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
            select_parents(arg, self.parents),
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

        def eq(a, b):
            if isinstance(a, float):
                return naneq(a, b)
            elif isinstance(a, np.ndarray):
                return np.array_equal(a, b)
            elif isinstance(a, Tuple):
                return all(eq(ai, bi) for (ai, bi) in zip(a, b))
            return a == b

        at = dataclasses.astuple(self)
        bt = dataclasses.astuple(other)
        return all(eq(a, b) for (a, b) in zip(at, bt))

    def copy(self):
        """
        Return event copy.
        """
        # Cannot be implemented with copy.deepcopy, leads to recursion
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
        return self._fast_selection(self.status == 1)

    def final_state_charged(self):
        """
        Return filtered event with only charged final state particles.

        The final state is generator-specific, see :meth:`final_state`.
        Additionally, this selects only charged particles which can be
        seen by a tracking detector.
        """
        return self._fast_selection((self.status == 1) & (self.charge != 0))

    def without_parton_shower(self):
        """
        Return filtered event without parton shower.

        Pythia generates a parton shower and makes it part of the history.
        Looking at this is interesting for generator experts, but for most
        analyses, this is not needed and can be removed from the particle
        history.
        """
        mask = True
        apid = np.abs(self.pid)
        for pid in quarks_and_diquarks_and_gluons:
            mask &= apid != pid
        return self[mask]

    def _fast_selection(self, arg):
        # This selection is faster than __getitem__, because we skip
        # parent selection, which is just wasting time if we select only
        # final state particles.
        save = self.parents
        self.parents = None
        event = self[arg]
        self.parents = save
        return event

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
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.log((self.p_tot + self.pz) / self.pt)

    @property
    def y(self):
        """Return rapidity."""
        with np.errstate(divide="ignore", invalid="ignore"):
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

        model, version = self.generator

        if genevent is None:
            genevent = pyhepmc.GenEvent()
            genevent.run_info = pyhepmc.GenRunInfo()
            genevent.run_info.tools = [(model, version, "")]

        # We must apply some workarounds so that HepMC3 conversion and IO works
        # for all models. This should be revisited once the fundamental issues
        # with particle histories have been fixed.
        if model == "Pythia" and version.startswith("8"):
            # to get a valid GenEvent we must
            # 1) select only particles produced after the parton shower
            # 2) connect particles attached to a single beam particle (diffractive events)
            #    to the common interaction vertex (1, 2)
            # TODO check if this costs significant amount of time and speed it up if so
            mask = (self.status == 1) | (self.status == 2) | (self.status == 4)
            ev = self[mask]
            mask = (ev.parents[:, 0] == 1) | (ev.parents[:, 0] == 2)
            ev.parents[mask] = (1, 2)
        elif model in ("UrQMD", "PhoJet", "DPMJET-III"):
            # can only save final state until history is fixed
            warnings.warn(
                f"{model}-{version}: only final state particles "
                "available in HepMC3 event",
                RuntimeWarning,
            )
            ev = self.final_state()
        else:
            ev = self

        genevent.from_hepevt(
            event_number=ev.nevent,
            px=ev.px,
            py=ev.py,
            pz=ev.pz,
            en=ev.en,
            m=ev.m,
            pid=ev.pid,
            status=ev.status,
            parents=ev.parents,
            children=ev.children,
            vx=ev.vx,
            vy=ev.vy,
            vz=ev.vz,
            vt=ev.vt,
        )

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

    def __init__(self, generator: MCRun, kinematics: EventKinematics):
        """
        Parameters
        ----------
        generator : Model
            Generator instance.
        kinematics : EventKinematics
            Kinematics of the event.
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

        EventData.__init__(
            self,
            (generator.name, generator.version),
            kinematics,
            int(getattr(evt, self._nevhep)),
            self._get_impact_parameter(),
            self._get_n_wounded(),
            getattr(evt, self._idhep)[sel],
            getattr(evt, self._isthep)[sel],
            self._charge_init(npart),
            *phep,
            *vhep,
            parents,
            children,
        )

    @abstractmethod
    def _charge_init(self, npart: int):
        # override this in derived to get charge info
        ...

    def _get_impact_parameter(self):
        # override this in derived
        return np.nan

    def _get_n_wounded(self):
        # override this in derived
        return (0, 0)

    # MCEvent is not pickleable, but EventData is. For convenience, we
    # make it so that MCEvent can be saved and is restored as EventData.
    def __getstate__(self):
        # Save only EventData sub-state. We actually need to make a copy here,
        # a view is not sufficient.
        copied_event = super().copy()
        return dataclasses.asdict(copied_event)

    def __getnewargs__(self):
        # upon unpickling, create EventData object instead of MCEvent object
        return (EventData,)


@dataclasses.dataclass
class RMMARDState:
    _c_number: np.ndarray
    _u_array: np.ndarray
    _u_i: np.ndarray
    _u_j: np.ndarray
    _seed: np.ndarray
    _counter: np.ndarray
    _big_counter: np.ndarray
    _sequence_number: np.ndarray
    _composite_target_state = None
    # only needed for Sibyll
    _gasdev_iset: np.ndarray = None

    def __init__(self, generator):
        data = generator._lib.crranma4

        self._c_number = data.c
        self._u_array = data.u
        self._u_i = data.i97
        self._u_j = data.j97
        self._seed = data.ijkl
        self._counter = data.ntot
        self._big_counter = data.ntot2
        self._sequence_number = data.jseq

        self._composite_target_state = generator._composite_target_rng.__getstate__()

        if generator.name == "SIBYLL":
            self._gasdev_iset = generator._lib.rndmgas.iset

    def _restore(self, generator):
        data = generator._lib.crranma4

        data.c = self._c_number
        data.u = self._u_array
        data.i97 = self._u_i
        data.j97 = self._u_j
        data.ijkl = self._seed
        data.ntot = self._counter
        data.ntot2 = self._big_counter
        data.jseq = self._sequence_number

        generator._composite_target_rng.__setstate__(self._composite_target_state)

        if generator.name == "SIBYLL":
            generator._lib.rndmgas.iset = self._gasdev_iset

    def __eq__(self, other):
        a = dataclasses.astuple(self)
        b = dataclasses.astuple(other)
        return all(np.all(ai == bi) for (ai, bi) in zip(a, b))

    def copy(self):
        """
        Return generator copy.
        """
        return copy.deepcopy(self)  # this uses setstate, getstate

    # HD: These should not be public, since the random number state
    # is an implementation detail. I am going to leave it, since this
    # whole class will become obsolete when we switch to the numpy PRNG.

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
    _once_called = False
    _alive_instances = set()
    _stable = set()

    # defaults for many models (override in Derived if needed)
    _projectiles = standard_projectiles
    _targets = Nuclei()
    _ecm_min = 10 * GeV

    nevents = 0  # number of generated events so far

    def __init__(self, seed, *args):
        import importlib
        from random import randint

        if self.pyname in self._alive_instances:
            warnings.warn(
                f"A previous instance of {self.pyname} is still alive. "
                "You cannot use two instances in parallel. "
                "Please delete the old one first before creating a new one. "
                "You can ignore this warning if the previous instance is already "
                "out of scope, Python does not always destroy old instances immediately.",
                RuntimeWarning,
                stacklevel=3,
            )

        self._alive_instances.add(self.pyname)

        if seed is None:
            self._seed = randint(1, 10000000)
        elif isinstance(seed, int):
            self._seed = seed
        else:
            raise ValueError(f"Invalid seed {seed}")

        # TODO use single PRNG for everything
        self._composite_target_rng = np.random.default_rng(self._seed)

        if not self._once_called:
            self._once_called = True
            assert hasattr(self, "_name")
            assert hasattr(self, "_version")
            assert hasattr(self, "_library_name")
            assert hasattr(self, "_event_class")
            assert hasattr(self, "_frame")
            self._lib = importlib.import_module(f"impy.models.{self._library_name}")
            if hasattr(self._lib, "init_rmmard"):
                self._lib.init_rmmard(self._seed)
            # Run internal model initialization code
            self._once(*args)
        else:
            if hasattr(self._lib, "init_rmmard"):
                self._lib.init_rmmard(self._seed)

        # Set standard long lived particles as stable
        for pid in long_lived:
            self.set_stable(pid)

    def __del__(self):
        if self.pyname in self._alive_instances:
            self._alive_instances.remove(self.pyname)

    def __call__(self, kin, nevents):
        """Generator function (in python sence)
        which launches the underlying event generator
        and returns its the result (event) as MCEvent object
        """
        nretries = 0
        for nev in self._composite_plan(kin, nevents):
            while nev > 0:
                if self._generate():
                    nretries = 0
                    self.nevents += 1
                    nev -= 1
                    event = self._event_class(self, kin)
                    # boost into frame requested by user
                    kin.apply_boost(event, self._frame)
                    yield event
                    continue
                nretries += 1
                if nretries % 50 == 0:
                    warnings.warn(
                        f"Event was rejected {nretries} times in a row. The generator "
                        "may be misconfigured or used outside of its valid range",
                        RuntimeWarning,
                    )
                if nretries > 1000:
                    raise RuntimeError("More than 1000 retries, aborting")

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
    def version(cls):
        """Event generator version"""
        return cls._version

    @classproperty
    def projectiles(cls):
        """Supported projectiles (positive PDGIDs only, c.c. implied)."""
        return cls._projectiles

    @classproperty
    def targets(cls):
        """Supported targets (positive PDGIDs only, c.c. implied)."""
        return cls._targets

    @classmethod
    def _validate_kinematics(cls, kin):
        if abs(kin.p1) not in cls._projectiles:
            raise ValueError(
                f"projectile {pdg2name(kin.p1)}[{int(kin.p1)}] is not allowed, "
                f"see {cls.pyname}.projectiles"
            )
        if abs(kin.p2) not in cls._targets:
            raise ValueError(
                f"target {pdg2name(kin.p2)}[{int(kin.p2)}] is not among allowed, "
                f"see {cls.pyname}.targets"
            )
        if kin.ecm < cls._ecm_min:
            raise ValueError(
                f"center-of-mass energy {kin.ecm/GeV} GeV < "
                f"minimum energy {cls._ecm_min/GeV} GeV"
            )

    def _composite_plan(self, kin, nevents):
        self._validate_kinematics(kin)
        if isinstance(kin.p2, CompositeTarget):
            nevents = self._composite_target_rng.multinomial(nevents, kin.p2.fractions)
            components = kin.p2.components
            # as a workaround for DPMJet, we generate heaviest elements first
            pairs = sorted(zip(components, nevents), key=lambda p: p[0].A, reverse=True)
            for c, nev in pairs:
                kin.p2 = c
                self._set_kinematics(kin)
                yield nev
        else:
            self._set_kinematics(kin)
            yield nevents

    @property
    def random_state(self):
        return RMMARDState(self)

    @random_state.setter
    def random_state(self, rng_state):
        rng_state._restore(self)

    def get_stable(self):
        """Return set of stable particles."""
        # make a copy to prevent modification of internal state
        return set(self._stable)

    def set_stable(self, particle, stable=True):
        """Prevent decay of an unstable particle.

        Parameters
        ----------
        particle : str or int
            Name or PDG ID of the particle.
        stable : bool, optional
            If true, particle is not decayed by the generator.

        Notes
        -----
        Some generators (e.g. the QGSJet models) do not provide
        a full particle history and do not allow one to set
        certain resonances as stable.
        """
        pid = process_particle(particle)
        p = Particle.from_pdgid(pid)
        if p.ctau is None or p.ctau == np.inf:
            raise ValueError(f"{pdg2name(pid)} cannot decay")
        if abs(pid) == 311:
            self.set_stable(130, stable)
            self.set_stable(310, stable)
            return

        if stable:
            self._stable.add(pid)
        elif pid in self._stable:
            self._stable.remove(pid)
        self._set_stable(pid, stable)

    def maydecay(self, particle):
        """Decay particle in event record.

        Equivalent to `self.set_stable(particle, stable=False)`

        Parameters
        ----------
        particle : str or int
            Name or PDG ID of the particle.
        """
        self.set_stable(particle, False)

    def cross_section(self, kin, **kwargs):
        """Cross sections according to current setup.

        Parameters
        ----------
        kin : EventKinematics
            Calculate cross-section for EventKinematics.
        kwargs :
            Further arguments passed to the model implementation.
        """
        self._validate_kinematics(kin)
        if isinstance(kin.p2, CompositeTarget):
            cross_section = CrossSectionData(0, 0, 0, 0, 0, 0, 0)
            components = kin.p2.components
            fractions = kin.p2.fractions
            for component, fraction in zip(components, fractions):
                kin.p2 = component
                cs = self._cross_section(kin, **kwargs)
                for i, val in enumerate(dataclasses.astuple(cs)):
                    cross_section[i] += fraction * val
            return cross_section
        else:
            return self._cross_section(kin, **kwargs)

    @abstractmethod
    def _once(self, *args):
        # This has to be implemented in the derived concrete class.
        pass

    @abstractmethod
    def _set_kinematics(self, kin):
        # Set new combination of energy, momentum, projectile
        # and target combination in the underlying model.

        # Either, this method defines some derived variables
        # that _generate() can use to generate new events
        # without additional arguments, or, it can also set
        # internal variables of the model. In both cases the
        # important thing is that _generate() remains argument-free.
        pass

    @abstractmethod
    def _set_stable(self, pid, stable):
        # This has to be implemented in the derived concrete class.
        pass

    @abstractmethod
    def _cross_section(self, kin):
        # This has to be implemented in the derived concrete class.
        pass

    @abstractmethod
    def _generate(self):
        # This has to be implemented in the derived concrete class.
        # Returns True if event was successfully generated and False otherwise.
        pass
