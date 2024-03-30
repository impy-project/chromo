"""The :mod:`common`module contains the classes that expose the
interaction model interface to the front-end and/or the user.

Classes derived from :class:`MCEvent` in (for example :mod:`chromo.models.sibyll`)
cast the data from the event generators' particle stacks to numpy arrays.
The basic variables are sufficient to compute all derived attributes,
such as the rapidity :func:`MCEvent.y` or the laboratory momentum fraction
:func:`MCEvent.xlab`.
"""

import copy
import dataclasses
import importlib
import warnings
from abc import ABC, abstractmethod
from contextlib import contextmanager
from typing import Optional, Tuple

import numpy as np
from packaging.version import parse as parse_version
from particle import Particle

from chromo.constants import (
    GeV,
    long_lived,
    quarks_and_diquarks_and_gluons,
    standard_projectiles,
)
from chromo.decay_handler import Pythia8DecayHandler
from chromo.kinematics import CompositeTarget, EventKinematicsBase
from chromo.util import (
    Nuclei,
    classproperty,
    naneq,
    pdg2name,
    select_long_lived,
    select_mothers,
    unique_sorted_pids,
)

all_unstable_pids = select_long_lived()


# Do we need EventData.n_spectators in addition to EventData.n_wounded?
# n_spectators can be computed from n_wounded as number_of_nucleons - n_wounded
# If we want this, it should be computed dynamically via a property.
@dataclasses.dataclass
class CrossSectionData:
    """Information of cross sections returned by the generator.

    The class is driven by the amount of detail that Pythia-8 provides.
    Most other generators do not fill all attributes. When information is missing,
    the attributes contain NaN. The fields are intentionally redundant, since
    some generators may provide only the total cross section or only the inelastic
    cross section.

    All cross sections are in millibarn.

    Attributes
    ----------
    total : float
        Total cross section (elastic + inelastic).
    inelastic : float
        Inelastic cross section. Includes diffractive cross sections.
    elastic : float
        Cross section for pure elastic scattering of incoming particle
        without any new particle generation.
    prod : float
        Particle production cross section, defined as total minus elastic.
    quasielastic : float
        Total quasielastic cross section (includes elastic).
    coherent: float
        (elastic with respect to the projectile) cross section
    diffractive_xb : float
        Single diffractive cross section. Particle 2 remains intact,
        particles are produced around Particle 1.
    diffractive_ax : float
        Single diffractive cross section. Particle 1 remains intact,
        particles are produced around Particle 2.
    diffractive_xx : float
        Double diffractive cross section. Central rapidity gap, but
        paricles are producted around both beam particles.
    diffractive_axb : float
        Central diffractive cross section (e.g.
        pomeron-pomeron interaction.)
    diffractive_sum : float
        Sum of diffractive cross sections.
    b_elastic : float
        Slope of elastic cross section in mb/GeV^2.
    """

    total: float = np.nan
    inelastic: float = np.nan
    elastic: float = np.nan
    prod: float = np.nan
    quasielastic: float = np.nan
    coherent: float = np.nan
    diffractive_xb: float = np.nan
    diffractive_ax: float = np.nan
    diffractive_xx: float = np.nan
    diffractive_axb: float = np.nan
    diffractive_sum: float = np.nan
    b_elastic: float = np.nan

    @property
    def non_diffractive(self):
        return self.inelastic - self.diffractive

    @property
    def diffractive(self):
        if float(self.diffractive_sum) > 0.0:
            return self.diffractive_sum
        return np.nansum(
            (
                self.diffractive_xb,
                self.diffractive_ax,
                self.diffractive_xx,
                self.diffractive_axb,
            )
        )

    def __eq__(self, other, rtol=1e-3):
        at = dataclasses.astuple(self)
        bt = dataclasses.astuple(other)
        return all(naneq(a, b, rtol) for (a, b) in zip(at, bt))

    def __ne__(self, other):
        return not self == other

    def _mul_radd(self, factor, other):
        for field in dataclasses.fields(self):
            setattr(
                self,
                field.name,
                getattr(self, field.name) + factor * getattr(other, field.name),
            )


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
    kin: EventKinematicsBase
        Info about initial state.
    nevent: int
        Which event in the sequence.
    impact_parameter: float
        Impact parameter for nuclear collisions in mm.
    n_wounded: (int, int)
        Number of wounded nucleons on sides A and B.
    production_cross_section: float
        Production cross section in mb; inelastic for photon-/hadron-hadron
        and production cross section for hadron-/nucleus-nucleus collisions.
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
    mothers: 2D array of int or None
        This array has shape (N, 2) for N particles. The two numbers
        are a range of indices pointing to parent particle(s). This
        array is not available for all generators and it is not
        available for filtered events. In those cases, mothers is None.
    daughters: 2D array of int or None
        Same as mothers.
    """

    generator: Tuple[str, str]
    kin: EventKinematicsBase
    nevent: int
    impact_parameter: float
    n_wounded: Tuple[int, int]
    production_cross_section: float
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
    mothers: Optional[np.ndarray]
    daughters: Optional[np.ndarray]

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
            self.production_cross_section,
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
            self.mothers.copy() if self.mothers is not None else None,
            self.daughters.copy() if self.daughters is not None else None,
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

        The final state contains all particles which do not have children.

        Unless the generator configuration is modified (see exceptions below),
        this means that the final state consists of prompt long-lived particles.

        Long-lived particles have live-times > 30 ps. Prompt particles either
        originate directly from primary interaction or have only parents with
        life-times < 30 ps. This is the ALICE definition of "promptness",
        see ALICE-PUBLIC-2017-005, which is identical to the definition by the
        "LHC Physics Center at CERN Minimum Bias and Underlying Event working group".

        Exceptions: Some generators deviate from this scheme. SIBYLL-2.1 does not
        produce Omega- and its antiparticle, so the final state never contains them.
        The QGSJet family does not produce Omega-, Xi-, Xi0, Sigma-, Sigma+ and their
        antiparticles.
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

    def _select(self, arg, update_mothers):
        # This selection is faster than __getitem__, because we skip
        # parent selection, which is just wasting time if we select only
        # final state particles.
        return EventData(
            self.generator,
            self.kin,
            self.nevent,
            self.impact_parameter,
            self.n_wounded,
            self.production_cross_section,
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
            select_mothers(arg, self.mothers) if update_mothers else None,
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

    @property
    def fw(self):
        """Quantity needed for invariant cross section histograms."""
        return self.en / self.kin.pcm

    def _prepare_for_hepmc(self):
        """
        Override this method in classes that need to modify event
        history for compatibility with HepMC.
        """
        return self

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

        if parse_version(pyhepmc.__version__) < parse_version("2.13.2"):
            raise RuntimeError(
                f"current pyhepmc version is {pyhepmc.__version__} < 2.13.2"
                f"\nPlease `pip install pyhepmc==2.13.2` or later version",
            )

        model, version = self.generator

        if genevent is None:
            genevent = pyhepmc.GenEvent()
            genevent.run_info = pyhepmc.GenRunInfo()
            genevent.run_info.tools = [(model, version, "")]

        ev = self._prepare_for_hepmc()
        genevent.from_hepevt(
            event_number=ev.nevent,
            px=ev.px,
            py=ev.py,
            pz=ev.pz,
            en=ev.en,
            m=ev.m,
            pid=ev.pid,
            status=ev.status,
            parents=ev.mothers if ev.mothers is not None else None,
            children=ev.daughters if ev.daughters is not None else None,
            vx=ev.vx,
            vy=ev.vy,
            vz=ev.vz,
            vt=ev.vt,
            fortran=False,
        )
        # Convert cross section from mb to pb
        genevent.cross_section = pyhepmc.GenCrossSection()
        genevent.cross_section.set_cross_section(
            self.production_cross_section * 1e9, 0, -1, -1
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

        self._generator_frame = generator._frame
        EventData.__init__(
            self,
            (generator.name, generator.version),
            generator.kinematics,
            int(getattr(evt, self._nevhep)),
            self._get_impact_parameter(),
            self._get_n_wounded(),
            generator._inel_or_prod_cross_section,
            getattr(evt, self._idhep)[sel],
            getattr(evt, self._isthep)[sel],
            self._charge_init(npart),
            *phep,
            *vhep,
            mothers=getattr(evt, self._jmohep).T[sel] if self._jmohep else None,
            daughters=getattr(evt, self._jdahep).T[sel] if self._jdahep else None,
        )

        if generator._restore_beam_and_history:
            self._history_zero_indexing()
            self._repair_initial_beam()

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

    def _repair_initial_beam(self):
        raise NotImplementedError("The method must be implemented in derived class")

    def _history_zero_indexing(self):
        self.mothers = self.mothers - 1
        self.daughters = self.daughters - 1

    def _prepend_initial_beam(self):
        beam = self.kin._get_beam_data(self._generator_frame)
        for field, beam_field in beam.items():
            event_field = getattr(self, field)
            if event_field is None:
                continue
            if field in ["mothers", "daughters"]:
                res = np.concatenate((beam_field, event_field + 2))
            else:
                res = np.concatenate((beam_field, event_field))
            setattr(self, field, res)


# =========================================================================
# MCRun
# =========================================================================
class MCRun(ABC):
    #: Prevent creating multiple classes within same python scope
    _is_initialized = []
    _restartable = False
    _projectiles = standard_projectiles
    _targets = Nuclei()
    _ecm_min = 10 * GeV  # default for many models
    # Corresponds to current cross section in mb, updated when kinematics is set
    _inel_or_prod_cross_section = None
    _restore_beam_and_history = True
    nevents = 0  # number of generated events so far
    _unstable_pids = set(all_unstable_pids)
    _final_state_particles = []
    _decay_handler = None  # Pythia8DecayHandler instance if activated

    def __init__(self, seed):
        if not self._restartable:
            self._abort_if_already_initialized()

        assert hasattr(self, "_name")
        assert hasattr(self, "_version")
        assert hasattr(self, "_library_name")
        assert hasattr(self, "_event_class")
        assert hasattr(self, "_frame")
        try:
            self._lib = importlib.import_module(f"chromo.models.{self._library_name}")
        except ModuleNotFoundError:
            self._lib = importlib.import_module(f"{self._library_name}")

        self._rng = np.random.default_rng(seed)
        if hasattr(self._lib, "npy"):
            self._lib.npy.bitgen = self._rng.bit_generator.ctypes.bit_generator.value

    def __call__(self, nevents):
        """Generator function (in python sence)
        which launches the underlying event generator
        and returns the event as MCEvent object
        """
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
                    self._validate_decay(event)
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
        # This is using private interface to get the seed. This may be brittle.
        # I did not find a way to get the seed using the public API.
        return self._rng._bit_generator._seed_seq.entropy

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

    def _check_kinematics(self, kin):
        """Check if kinematics are allowed for this generator."""

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

    @property
    def kinematics(self):
        return self._kinematics

    @kinematics.setter
    def kinematics(self, kin):
        self._check_kinematics(kin)
        self._kinematics = kin
        self._set_kinematics(kin)

        if (kin.p1.is_nucleus and kin.p1.A > 1) or (kin.p2.is_nucleus and kin.p2.A > 1):
            self._inel_or_prod_cross_section = self.cross_section().prod
        else:
            self._inel_or_prod_cross_section = self.cross_section().inelastic

    def cross_section(self, kin=None, max_info=False):
        """Cross sections according to current setup.

        Parameters
        ----------
        kin : EventKinematics, optional
            If provided, calculate cross-section for EventKinematics.
            Otherwise return values for current setup.
        max_info : bool, optional
            Return full maximal information about interaction cross sections for
            nucleus-nucleus case. Slow - uses Monte Carlo integration. Number of
            trials is controled separately with `generator.n_trials` attribute.
        """
        with self._temporary_kinematics(kin):
            kin2 = self.kinematics
            if isinstance(kin2.p2, CompositeTarget):
                cross_section = CrossSectionData(0, 0, 0, 0, 0, 0, 0)
                kin3 = copy.copy(kin2)
                for component, fraction in zip(kin2.p2.components, kin2.p2.fractions):
                    kin3.p2 = component
                    # this calls cross_section recursively, which is fine
                    cross_section._mul_radd(
                        fraction, self._cross_section(kin3, max_info=max_info)
                    )
                return cross_section
            else:
                return self._cross_section(kin, max_info=max_info)

    @abstractmethod
    def _cross_section(self, kin):
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

    def _validate_decay(self, event):
        """Checks whether all unstable particles are final state particles.
        If any unstable particles are not yet decayed, it attempts to decay them
        using the decay_handler.
        """
        final_pids = event.pid[event.status == 1]
        may_decay = np.isin(final_pids, all_unstable_pids)

        if not np.all(np.isin(final_pids[may_decay], self._final_state_particles)):
            if self._decay_handler:
                self._decay_handler(event)
            else:
                final_pids_dec = final_pids[may_decay]
                not_decayed = np.logical_not(
                    np.isin(final_pids_dec, self._final_state_particles)
                )
                not_decayed_pids = set(final_pids_dec[not_decayed])
                warnings.warn(
                    f"{self.pyname}: {not_decayed_pids} haven't been decayed. "
                    "Consider to use generator._activate_decay_handler(on=True)",
                    RuntimeWarning,
                )

    @property
    def final_state_particles(self):
        """Returns a sorted list of particles that can decay
        but are considered stable by the event generator."""
        return tuple(unique_sorted_pids(self._final_state_particles))

    @final_state_particles.setter
    def final_state_particles(self, list_of_pdgids):
        """
        Sets particles with PDG IDs provided in `pdgid` list as stable particles.
        All other unstable particles will decay. Anti-particles are synchronized.

        Stable particles in `pdgid` with tau0 = inf are ignored.

        Args:
            pdgid (list): A list of PDG IDs for particles that should be considered
                        stable (present in the final state).

        Example:
            To configure an `QGSJetII04` event generator to treat charged pions
            (PDG ID 211 and -211) and muons (PDG ID 13 and -13) as stable
            particles in the final state:

            >>> evt_kin = chromo.kinematics.FixedTarget(100, "p", "p")
            >>> generator = chromo.models.QGSJetII04(evt_kin)
            >>> generator.final_state_particles = [211, -211, 13, -13]

            If you need to set particles as stable for those with a lifetime
            greater than `tau_stable`:
            >>> generator.final_state_particles = (chromo.util
                                                  .select_long_lived(tau_stable))
        """
        self._set_final_state_particles(list_of_pdgids)

    def _set_final_state_particles(self, list_of_pdgids=long_lived):
        """By default defines particles as stable
        for the default 'tau_stable' value in the config."""

        fsparticles = np.unique(list_of_pdgids)
        is_unstable = np.isin(fsparticles, list(self._unstable_pids))
        fsparticles = fsparticles[is_unstable]
        # Clean up by setting all unstable particles as unstable
        for pid in self._unstable_pids:
            self.set_stable(pid, False, update_decay_handler=False)

        for pid in fsparticles:
            self.set_stable(pid, True, update_decay_handler=False)

        self._sync_decay_handler()

    def set_stable(self, pdgid, stable=True, update_decay_handler=True):
        """Prevent decay of unstable particles.

        Anti-particles are synchronized.

        Args:
            pdgid (int)        : PDG ID of the particle
            stable (bool)      : If `False`, particle is allowed to decay
        """
        p = Particle.from_pdgid(pdgid)
        assert pdgid in self._unstable_pids, f"{pdg2name(pdgid)} unknown or stable"
        ap = p.invert() if p.invert() != p else False
        if p.ctau is None or p.ctau == np.inf:
            raise ValueError(f"{pdg2name(pdgid)} cannot decay")

        if abs(pdgid) == 311:
            pdgid_list = [130, 310]
        elif ap:
            pdgid_list = [pdgid, ap.pdgid]
        else:
            pdgid_list = [pdgid]

        for pdgid in pdgid_list:
            self._set_stable(pdgid, stable)

        if stable:
            self._final_state_particles = np.unique(
                np.append(self._final_state_particles, pdgid_list).astype(np.int64)
            )
        else:
            if len(self._final_state_particles) > 0:
                remove = np.isin(self._final_state_particles, pdgid_list)
                self._final_state_particles = self._final_state_particles[~remove]

        if update_decay_handler:
            self._sync_decay_handler()

    def _sync_decay_handler(self):
        # Synchronize decay handler
        if self._decay_handler:
            self._decay_handler.set_stable(self._final_state_particles)
            assert np.isin(
                self._final_state_particles, self._decay_handler.all_stable_pids
            ).all(), "Decay handler and generator are out of sync"

    def set_unstable(self, pdgid):
        """Convenience funtion for `self.set_stable(..., stable=False)`

        Args:
            pdgid(int)         : PDG ID of the particle
        """
        self.set_stable(pdgid, False)

    @abstractmethod
    def _set_stable(self, pdgid, stable):
        pass

    def _activate_decay_handler(self, on=True, *, seed=None):
        """
        Activates the Pythia8 decay handler for any of the generators
        except Pythia8 itself.

        This function is responsible for activating the decay handler which
        ensures that particles, which are set to be unstable actually decay
        consistently. This feature mainly fixes lacking decay functions in
        models like QGSJet. There is some notable but not dramatic performance
        impact.

        The function is private since the stability and the interface is not
        guaranteed to last, and it doesn't work on Windows due to compilation
        issues of Pythia8 on that OS.

        Args:
            on (bool)       : If `True`, the decay handler is activated or destroyed

        Returns:
            None
        """
        if (not on) or (self.pyname == "Pythia8"):
            if self._decay_handler:
                del self._decay_handler
            self._decay_handler = None
            return

        if not self._decay_handler:
            try:
                self._decay_handler = Pythia8DecayHandler(
                    self._final_state_particles, seed=seed
                )
            except ModuleNotFoundError as ex:
                import warnings

                warnings.warn(
                    f"Pythia8DecayHandler not available:\n{ex}\n"
                    "Some particles may not decay",
                    RuntimeWarning,
                )
                self._decay_handler = None
                return

    @contextmanager
    def _temporary_kinematics(self, kin):
        if kin is None:
            yield
        else:
            prev = copy.copy(self.kinematics)
            self.kinematics = kin
            yield
            self.kinematics = prev
