""" This module handles transformations between Lorentz frames and
different inputs required by the low-level event generator interfaces.
"""
import numpy as np
from impy.util import (
    TaggedFloat,
    energy2momentum,
    momentum2energy,
    elab2ecm,
    ecm2elab,
    mass,
    is_real_nucleus,
    process_particle,
)
from impy.constants import nucleon_mass, MeV, GeV, TeV, PeV, EeV
from impy.util import CompositeTarget, EventFrame
from particle import PDGID
import dataclasses
from typing import Union, Tuple


__all__ = (
    "EventFrame",
    "CompositeTarget",
    "MeV",
    "GeV",
    "TeV",
    "PeV",
    "EeV",
    "EventKinematics",
    "CenterOfMass",
    "FixedTarget",
    "TotalEnergy",
    "KinEnergy",
    "Momentum",
)


@dataclasses.dataclass
class EventKinematics:
    """Handles kinematic variables and conversions between reference frames.

    There are different ways to specify a particle collision. For instance
    the projectile and target momenta can be specified in the target rest frame,
    the so called 'laboratory' frame, or the nucleon-nucleon center-of-mass frame
    where the modulus of the nucleon momenta is the same but the direction
    inverted. Each event generator expects its arguments to be given in one
    or the other frame. This class allows the generator to pick itself the correct
    frame, while the user can specify the kinematics in the preferred form.

    Parameters
    ----------
        particle1: str or int or (int, int)
            Particle name, PDG ID, or nucleus mass & charge (A, Z) of projectile.
        particle2: str or int or (int, int) or CompositeTarget
            Particle name, PDG ID, nucleus mass & charge (A, Z), or CompositeTarget
            of the target
        ecm : float, optional
            Center-of-mass energy :math:`\\sqrt{s}`.
        plab : float, optional
            Projectile momentum in lab frame. If the projectile is a nucleus, it is
            the momentum per nucleon.
        elab : float, optional
            Projectile energy in lab frame. If the projectile is a nucleus, it is
            the energy per nucleon.
        ekin : float, optional
            Projectile kinetic energy in lab frame. If the projectile is a nucleus,
            it is the kinetic energy per nucleon.
        beam : tuple of two floats
            Specification as tuple of two momenta. If the projectile or target are
            nuclei, it is the momentum per nucleon.
    """

    frame: EventFrame
    p1: Union[PDGID, Tuple[int, int]]
    p2: Union[PDGID, Tuple[int, int], CompositeTarget]
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

        self.p1 = process_particle(particle1)
        self.p2 = process_particle(particle2)

        if isinstance(self.p1, CompositeTarget):
            raise ValueError("Only 2nd particle can be CompositeTarget")

        p2_is_composite = isinstance(self.p2, CompositeTarget)

        m1 = nucleon_mass if is_real_nucleus(self.p1) else mass(self.p1)
        m2 = nucleon_mass if is_real_nucleus(self.p2) else mass(self.p2)

        self.beams = (np.zeros(4), np.zeros(4))

        # Input specification in center-of-mass frame
        if ecm is not None:
            self.frame = frame or EventFrame.CENTER_OF_MASS
            self.ecm = ecm
            self.elab = ecm2elab(ecm, m1, m2)
            self.plab = energy2momentum(self.elab, m1)
        # Input specification as 4-vectors
        elif beam is not None:
            if p2_is_composite:
                raise ValueError("beam cannot be used with CompositeTarget")
            self.frame = frame or EventFrame.GENERIC
            p1, p2 = beam
            self.beams[0][2] = p1
            self.beams[1][2] = p2
            self.beams[0][3] = momentum2energy(p1, m1)
            self.beams[1][3] = momentum2energy(p2, m2)
            s = np.sum(self.beams, axis=0)
            # We compute ecm with energy2momentum. It is not really energy to momentum,
            # but energy2momentum(x, y) computes x^2 - y^2, which is what we need. Here,
            # I use that px and py are always zero, if we ever change this, many formulas
            # have to change in this class, like all the boosts
            self.ecm = energy2momentum(s[3], s[2])
            self.elab = ecm2elab(self.ecm, m1, m2)
            self.plab = energy2momentum(self.elab, m1)
        # Input specification in lab frame
        elif elab is not None:
            if not (elab > m1):
                raise ValueError("projectile energy > projectile mass required")
            self.frame = frame or EventFrame.FIXED_TARGET
            self.elab = elab
            self.plab = energy2momentum(self.elab, m1)
            self.ecm = elab2ecm(self.elab, m1, m2)
        elif ekin is not None:
            self.frame = frame or EventFrame.FIXED_TARGET
            self.elab = ekin + m1
            self.plab = energy2momentum(self.elab, m1)
            self.ecm = elab2ecm(self.elab, m1, m2)
        elif plab is not None:
            self.frame = frame or EventFrame.FIXED_TARGET
            self.plab = plab
            self.elab = momentum2energy(self.plab, m1)
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
        # set energies
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

    def __repr__(self):
        p = f"{self.p1:d}, {self.p2:d}"
        if self.frame == EventFrame.CENTER_OF_MASS:
            s = f"ecm={self.ecm}"
            return f"CenterOfMass({self.ecm}, {p})"
        elif self.frame == EventFrame.FIXED_TARGET:
            return f"FixedTarget({self.elab}, {p})"
        s = f"beam=({self.beams[0][2]}, {self.beams[1][2]}), frame={self.frame}"
        return f"EventKinematics({p}, {s})"


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
