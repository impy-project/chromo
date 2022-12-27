""" This module handles transformations between Lorentz frames and
different inputs required by the low-level event generator interfaces.


@Hans @Sonia: we need to come up with some sort general handling
of different inputs. Hans suggested to mimic a similar behavior as for
colors in matplotlib. That one can initialize with different arguments, like
'pp' 7 TeV would be internally translated to 2212, 2212 and to 4-vectors
[0.,0.,+-3.499999TeV, 3.5TeV]. This assumes that the input interpreted as
center of mass total energy (not momentum) AND that the final state is
defined in center-of-mass as well.

This was already the initial motivation but I have the impression that
the implementation is very "cooked up". We have to discuss this.

"""

import numpy as np
from impy.util import (
    TaggedFloat,
    AZ2pdg,
    is_AZ,
    energy2momentum,
    elab2ecm,
    mass,
    pdg2name,
    name2pdg,
)
from impy.constants import nucleon_mass
from particle import PDGID
import dataclasses
from typing import Union, Tuple
from enum import Enum


EventFrame = Enum("EventFrame", ["CENTER_OF_MASS", "FIXED_TARGET", "GENERIC"])


@dataclasses.dataclass(init=False)
class CompositeTarget:
    """Definition of composite targets made of multiple (atomic) nuclei.

    Examples of such composite targets are Air, CO_2, HCl, C_2H_60.
    """

    label: str
    components: Tuple[PDGID]
    fractions: np.ndarray

    # TODO: this is missing a docstring
    def __init__(self, components, label=""):
        if len(components) == 0:
            raise ValueError("components cannot be empty")
        fractions = np.empty(len(components))
        c = []
        for i, (particle, amount) in enumerate(components):
            fractions[i] = amount
            p = _normalize_particle(particle)
            if not p.is_nucleus:
                raise ValueError(f"component {particle} is not a nucleus")
            c.append(p)
        self.label = label
        self.components = tuple(c)
        self.fractions = fractions / np.sum(fractions)
        self.fractions.flags["WRITEABLE"] = False

    @property
    def Z(self):
        """Return maximum charge number."""
        # needed for compatibility with PDGID interface and for dpmjet initialization
        return max(p.Z for p in self.components)

    @property
    def A(self):
        """Return maximum number of nucleons."""
        # needed for compatibility with PDGID interface and for dpmjet initialization
        return max(p.A for p in self.components)

    @property
    def is_nucleus(self):
        return True

    def __int__(self):
        """Return PDGID for heaviest of elements."""
        return int(max((c.A, c) for c in self.components)[1])

    def average_mass(self):
        return sum(
            f * p.A * nucleon_mass for (f, p) in zip(self.fractions, self.components)
        )

    def __abs__(self):
        return abs(int(self))

    def __repr__(self):
        components = [
            (pdg2name(c), amount)
            for (c, amount) in zip(self.components, self.fractions)
        ]
        args = f"{components}"
        if self.label:
            args += f", label={self.label!r}"
        return f"CompositeTarget({args})"


def _normalize_particle(x):
    if isinstance(x, (PDGID, CompositeTarget)):
        return x
    if isinstance(x, int):
        return PDGID(x)
    if isinstance(x, str):
        try:
            return PDGID(name2pdg(x))
        except KeyError:
            raise ValueError(f"particle with name {x} not recognized")
    if is_AZ(x):
        return PDGID(AZ2pdg(*x))
    raise ValueError(f"{x} is not a valid particle specification")


@dataclasses.dataclass
class EventKinematics:
    """Handles kinematic variables and conversions between reference frames.

    There are different ways to specify a particle collision. For instance
    the projectile and target momenta can be specified in the target rest frame,
    the so called 'laboratory' frame, or the nucleon-nucleon center of mass frame
    where the modulus of the nucleon momenta is the same but the direction
    inverted. Each event generator expects its arguments to be given in one
    or the other frame. This class allows the generator to pick itself the correct
    frame, while the user can specify the kinematics in the preferred form.

    Args:
        (float) ecm      : :math:`\\sqrt{s}`, the center-of-mass energy
        (float) plab     : projectile momentum in lab frame
        (float) elab     : projectile energy in lab frame
        (float) ekin     : projectile kinetic energy in lab frame
        (tuple) beam     : Specification as tuple of 4-vectors (np.array)s
        (tuple) particle1: particle name, PDG ID, or nucleus mass & charge (A, Z)
                           of the projectile
        (tuple) particle2: particle name, PDG ID, or nucleus mass & charge (A, Z),
                           or CompositeTarget of the target

    """

    frame: EventFrame
    p1: PDGID
    p2: Union[PDGID, CompositeTarget]
    ecm: float  # for ions this is nucleon-nucleon collision system
    beams: Tuple[np.ndarray, np.ndarray]
    _gamma_cm: float
    _betagamma_cm: float

    def __init__(
        self,
        particle1,
        particle2,
        *,
        ecm=None,
        plab=None,
        elab=None,
        ekin=None,
        beam=None,
        frame=None,
    ):
        # Catch input errors

        if sum(x is not None for x in [ecm, plab, elab, ekin, beam]) != 1:
            raise ValueError(
                "Please provide only one of ecm/plab/elab/ekin/beam arguments"
            )

        if particle1 is None or particle2 is None:
            raise ValueError("particle1 and particle2 must be set")

        self.p1 = _normalize_particle(particle1)
        self.p2 = _normalize_particle(particle2)

        if isinstance(self.p1, CompositeTarget):
            raise ValueError("Only 2nd particle can be CompositeTarget")

        p2_is_composite = isinstance(self.p2, CompositeTarget)

        m1 = mass(self.p1)
        m2 = nucleon_mass if p2_is_composite or self.p2.is_nucleus else mass(self.p2)

        self.beams = (np.zeros(4), np.zeros(4))

        # Input specification in center of mass frame
        if ecm is not None:
            self.frame = EventFrame.CENTER_OF_MASS if frame is None else frame
            self.ecm = ecm
            self.elab = 0.5 * (ecm**2 - m1**2 - m2**2) / m2
            self.plab = energy2momentum(self.elab, m1)
        # Input specification as 4-vectors
        elif beam is not None:
            if p2_is_composite:
                raise ValueError("beam cannot be used with CompositeTarget")
            self.frame = EventFrame.GENERIC if frame is None else frame
            p1, p2 = beam
            self.beams[0][2] = p1
            self.beams[1][2] = p2
            self.beams[0][3] = np.sqrt(m1**2 + p1**2)
            self.beams[1][3] = np.sqrt(m2**2 + p2**2)
            s = np.sum(self.beams, axis=0)
            self.ecm = np.sqrt(s[3] ** 2 - np.sum(s[:3] ** 2))
            self.elab = 0.5 * (self.ecm**2 - m1**2 + m2**2) / m2
            self.plab = energy2momentum(self.elab, m1)
        # Input specification in lab frame
        elif elab is not None:
            if not (elab > m1):
                raise ValueError("projectile energy > projectile mass required")
            self.frame = EventFrame.FIXED_TARGET if frame is None else frame
            self.elab = elab
            self.plab = energy2momentum(self.elab, m1)
            self.ecm = np.sqrt(2.0 * self.elab * m2 + m2**2 + m1**2)
            # self.ecm = np.sqrt((self.elab + m2)**2 - self.plab**2)
        elif ekin is not None:
            self.frame = EventFrame.FIXED_TARGET if frame is None else frame
            self.elab = ekin + m1
            self.plab = energy2momentum(self.elab, m1)
            self.ecm = elab2ecm(self.elab, m1, m2)
        elif plab is not None:
            self.frame = EventFrame.FIXED_TARGET if frame is None else frame
            self.plab = plab
            self.elab = np.sqrt(self.plab**2 + m1**2)
            self.ecm = elab2ecm(self.elab, m1, m2)
        else:
            assert False  # this should never happen

        self._fill_beams(m1, m2)

        self._gamma_cm = (self.elab + m2) / self.ecm
        self._betagamma_cm = self.plab / self.ecm

    def apply_boost(self, event, generator_frame):
        if generator_frame == self.frame:
            return
        CMS = EventFrame.CENTER_OF_MASS
        FT = EventFrame.FIXED_TARGET
        if generator_frame == FT and self.frame == CMS:
            bg = -self._betagamma_cm
        elif generator_frame == CMS and self.frame == FT:
            bg = self._betagamma_cm
        else:
            raise NotImplementedError(
                f"Boosts from {generator_frame} to {self.frame} are not yet supported"
            )
        g = self._gamma_cm
        en = g * event.en + bg * event.pz
        pz = bg * event.en + g * event.pz
        event.en[:] = en
        event.pz[:] = pz

    def _fill_beams(self, m1, m2):
        if self.frame == EventFrame.GENERIC:
            return  # nothing to do
        if self.frame == EventFrame.CENTER_OF_MASS:
            s = self.ecm**2
            pcm = np.sqrt((s - (m1 + m2) ** 2) * (s - (m1 - m2) ** 2)) / (2 * self.ecm)
            self.beams[0][2] = pcm
            self.beams[1][2] = -pcm
        elif self.frame == EventFrame.FIXED_TARGET:
            self.beams[0][2] = self.plab
            self.beams[1][2] = 0
        # set energyies
        for b, m in zip(self.beams, (m1, m2)):
            b[3] = np.sqrt(m**2 + b[2] ** 2)

    def __eq__(self, other):
        at = dataclasses.astuple(self)
        bt = dataclasses.astuple(other)

        def eq(a, b):
            if isinstance(a, Tuple):
                return all(eq(ai, bi) for (ai, bi) in zip(a, b))
            if isinstance(a, np.ndarray):
                return np.array_equal(a, b)
            return a == b

        return all(eq(a, b) for (a, b) in zip(at, bt))


class CenterOfMass(EventKinematics):
    def __init__(self, ecm, particle1, particle2):
        super().__init__(ecm=ecm, particle1=particle1, particle2=particle2)


class TotalEnergy(TaggedFloat):
    pass


class KinEnergy(TaggedFloat):
    pass


class Momentum(TaggedFloat):
    pass


class FixedTarget(EventKinematics):
    def __init__(self, energy, particle1, particle2):
        if isinstance(energy, (TotalEnergy, int, float)):
            super().__init__(
                elab=float(energy), particle1=particle1, particle2=particle2
            )
        elif isinstance(energy, KinEnergy):
            super().__init__(
                ekin=float(energy), particle1=particle1, particle2=particle2
            )
        elif isinstance(energy, Momentum):
            super().__init__(
                plab=float(energy), particle1=particle1, particle2=particle2
            )
        else:
            raise ValueError(
                f"{energy!r} is neither a number nor one of "
                "TotalEnergy, KinEnergy, Momentum"
            )
