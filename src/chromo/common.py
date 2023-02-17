"""The :mod:`common`module contains the classes that expose the
interaction model interface to the front-end and/or the user.

Classes derived from :class:`MCEvent` in (for example :mod:`chromo.models.sibyll`)
cast the data from the event generators' particle stacks to numpy arrays.
The basic variables are sufficient to compute all derived attributes,
such as the rapidity :func:`MCEvent.y` or the laboratory momentum fraction
:func:`MCEvent.xlab`.
"""
from abc import ABC, abstractmethod
import numpy as np
from chromo.util import (
    classproperty,
    select_parents,
    naneq,
    pdg2name,
    Nuclei,
)
from chromo.constants import (
    quarks_and_diquarks_and_gluons,
    long_lived,
    standard_projectiles,
    GeV,
)
from chromo.kinematics import EventKinematics, CompositeTarget
import dataclasses
import copy
from typing import Tuple, Optional
from contextlib import contextmanager
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
        return self._select(arg, True)

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

    def __getstate__(self):
        t = [
            self.generator,
            self.kin.copy(),
            self.nevent,
            self.impact_parameter,
            self.n_wounded,
            self.pid.copy(),
            self.status.copy(),
            self.charge.copy(),
            self.px.copy(),
            self.py.copy(),
            self.pz.copy(),
            self.en.copy(),
            self.m.copy(),
            self.vx.copy(),
            self.vy.copy(),
            self.vz.copy(),
            self.vt.copy(),
            self.parents.copy(),
            self.children.copy() if self.children is not None else None,
        ]
        return t

    def __setstate__(self, state):
        for f, v in zip(dataclasses.fields(self), state):
            setattr(self, f.name, v)

    def copy(self):
        """
        Return event copy.
        """
        copies = self.__getstate__()
        return EventData(*copies)

    def final_state(self):
        """
        Return filtered event with only final state particles.

        The final state is generator-specific, but usually contains only
        long-lived particles.
        """
        return self._select(self.status == 1, False)

    def final_state_charged(self):
        """
        Return filtered event with only charged final state particles.

        The final state is generator-specific, see :meth:`final_state`.
        Additionally, this selects only charged particles which can be
        seen by a tracking detector.
        """
        return self._select((self.status == 1) & (self.charge != 0), False)

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

    def _select(self, arg, update_parents):
        # This selection is faster than __getitem__, because we skip
        # parent selection, which is just wasting time if we select only
        # final state particles.
        pid = self.pid[arg]
        return EventData(
            self.generator,
            self.kin,
            self.nevent,
            self.impact_parameter,
            self.n_wounded,
            pid,
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
            select_parents(arg, self.parents) if update_parents else None,
            None,
        )

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
        from chromo.util import EventFrame

        kin = self.kin
        if kin.frame == EventFrame.FIXED_TARGET:
            return self.en
        return kin._gamma_cm * self.en + kin._betagamma_cm * self.pz

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

        EventData.__init__(
            self,
            (generator.name, generator.version),
            generator.kinematics,
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
    def _charge_init(self, npart):
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
        # save only EventData sub-state
        return EventData.__getstate__(self)

    def __getnewargs__(self):
        # upon unpickling, create EventData object instead of MCEvent object
        return (EventData,)


# =========================================================================
# MCRun
# =========================================================================
class MCRun(ABC):
    #: Prevent creating multiple classes within same python scope
    _is_initialized = []
    _restartable = False
    _set_final_state_particles_called = False
    _projectiles = standard_projectiles
    _targets = Nuclei()
    _ecm_min = 10 * GeV  # default for many models
    nevents = 0  # number of generated events so far

    def __init__(self, seed):
        import importlib

        if not self._restartable:
            self._abort_if_already_initialized()

        assert hasattr(self, "_name")
        assert hasattr(self, "_version")
        assert hasattr(self, "_library_name")
        assert hasattr(self, "_event_class")
        assert hasattr(self, "_frame")
        self._lib = importlib.import_module(f"chromo.models.{self._library_name}")

        self._seed = seed
        self._rng = np.random.default_rng(seed)
        if hasattr(self._lib, "npy"):
            self._lib.npy.bitgen = self._rng.bit_generator.ctypes.bit_generator.value

    def __call__(self, nevents):
        """Generator function (in python sence)
        which launches the underlying event generator
        and returns its the result (event) as MCEvent object
        """
        assert self._set_final_state_particles_called
        nretries = 0
        for nev in self._composite_plan(nevents):
            while nev > 0:
                if self._generate():
                    nretries = 0
                    self.nevents += 1
                    nev -= 1
                    event = self._event_class(self)
                    # boost into frame requested by user
                    self.kinematics.apply_boost(event, self._frame)
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

    @abstractmethod
    def _generate(self):
        """The method to generate a new event.

        Returns
        -------
        bool
            True if event was successfully generated and False otherwise.
        """
        pass

    @abstractmethod
    def _set_kinematics(self, kin):
        # Set new combination of energy, momentum, projectile
        # and target combination for next event.

        # Either, this method defines some derived variables
        # that _generate_event() can use to generate new events
        # without additional arguments, or, it can also set
        # internal variables of the model. In both cases the
        # important thing is that _generate_event remains argument-free.

        # This call must not update self.kinematics, only the
        # generator-specific variables.
        pass

    def _composite_plan(self, nevents):
        kin = self.kinematics
        if isinstance(kin.p2, CompositeTarget):
            nevents = self._rng.multinomial(nevents, kin.p2.fractions)
            ek = copy.deepcopy(kin)
            for c, k in zip(kin.p2.components, nevents):
                ek.p2 = c
                with self._temporary_kinematics(ek):
                    yield k
        else:
            yield nevents

    @property
    def random_state(self):
        return self._rng.__getstate__()

    @random_state.setter
    def random_state(self, rng_state):
        self._rng.__setstate__(rng_state)

    @property
    def kinematics(self):
        return self._kinematics

    @kinematics.setter
    def kinematics(self, kin):
        if abs(kin.p1) not in self._projectiles:
            raise ValueError(
                f"projectile {pdg2name(kin.p1)}[{int(kin.p1)}] is not allowed, "
                f"see {self.pyname}.projectiles"
            )
        if abs(kin.p2) not in self._targets:
            raise ValueError(
                f"target {pdg2name(kin.p2)}[{int(kin.p2)}] is not among allowed, "
                f"see {self.pyname}.targets"
            )
        if kin.ecm < self._ecm_min:
            raise ValueError(
                f"center-of-mass energy {kin.ecm/GeV} GeV < "
                f"minimum energy {self._ecm_min/GeV} GeV"
            )
        self._kinematics = kin
        self._set_kinematics(kin)

    def set_stable(self, pdgid, stable=True):
        """Prevent decay of unstable particles

        Args:
            pdgid (int)        : PDG ID of the particle
            stable (bool)      : If `False`, particle is allowed to decay
        """
        p = Particle.from_pdgid(pdgid)
        if p.ctau is None or p.ctau == np.inf:
            raise ValueError(f"{pdg2name(pdgid)} cannot decay")
        if abs(pdgid) == 311:
            self._set_stable(130, stable)
            self._set_stable(310, stable)
        else:
            self._set_stable(pdgid, stable)

    def set_unstable(self, pdgid):
        """Convenience funtion for `self.set_stable(..., stable=False)`

        Args:
            pdgid(int)         : PDG ID of the particle
        """
        self.set_stable(pdgid, False)

    def cross_section(self, kin=None):
        """Cross sections according to current setup.

        Parameters
        ----------
        kin : EventKinematics, optional
            If provided, calculate cross-section for EventKinematics.
            Otherwise return values for current setup.
        """
        with self._temporary_kinematics(kin):
            kin2 = self.kinematics
            if isinstance(kin2.p2, CompositeTarget):
                cross_section = CrossSectionData(0, 0, 0, 0, 0, 0, 0)
                kin3 = copy.copy(kin2)
                for component, fraction in zip(kin2.p2.components, kin2.p2.fractions):
                    kin3.p2 = component
                    # this calls cross_section recursively, which is fine
                    cs = self.cross_section(kin3)
                    for field, value in dataclasses.asdict(cs).items():
                        val = getattr(cross_section, field) + fraction * value
                        setattr(cross_section, field, val)
                return cross_section
            else:
                return self._cross_section(kin)

    @abstractmethod
    def _cross_section(self, kin):
        pass

    @abstractmethod
    def _set_stable(self, pidid, stable):
        pass

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

    def _set_final_state_particles(self):
        """Defines particles as stable for the default 'tau_stable'
        value in the config."""

        for pdgid in long_lived:
            self._set_stable(pdgid, True)

        self._set_final_state_particles_called = True

    @contextmanager
    def _temporary_kinematics(self, kin):
        if kin is None:
            yield
        else:
            prev = copy.copy(self.kinematics)
            self.kinematics = kin
            yield
            self.kinematics = prev
