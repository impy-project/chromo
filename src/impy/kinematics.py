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
from impy import pdata, impy_config
from impy.util import info
import particle


def _is_AZ_tuple(arg):
    if isinstance(arg, tuple):
        if len(arg) != 2:
            raise ValueError(f"tuple (A, Z) should have len == 2, but given = {arg}")
        if not isinstance(arg[0], int):
            raise ValueError(
                "1st entry of (A, Z) should be 'int' type, but it is "
                f"{type(arg[0])} = {arg[0]}"
            )
        if not isinstance(arg[1], int):
            raise ValueError(
                "2nd entry of (A, Z) should be 'int' type, but it is "
                f"{type(arg[1])} = {arg[1]}"
            )
        return True
    else:
        return False


class _FromParticleName:
    all_pdgs = {p.name: int(p.pdgid) for p in particle.Particle.findall()}
    all_pdgs.update(
        {p.programmatic_name: int(p.pdgid) for p in particle.Particle.findall()}
    )
    all_pdgs.update(
        photon=22,
        Higgs=25,
        proton=2212,
        antiproton=-2212,
        neutron=2112,
        antineutron=-2112,
    )

    @staticmethod
    def _get_pdg(pname):
        try:
            return _FromParticleName.all_pdgs[pname]
        except KeyError:
            raise ValueError(f'Particle with name = "{pname}" is not found')

    @staticmethod
    def _get_AZ(arg):
        if isinstance(arg, str):
            pdg = particle.pdgid.PDGID(_FromParticleName._get_pdg(arg))
            if (pdg.A is None) or (pdg.Z is None):
                raise ValueError(f"'_get_AZ': no (A, Z) data for '{arg}'")
            else:
                return pdg.A, pdg.Z
        elif _is_AZ_tuple(arg):
            return arg
        else:
            raise ValueError(
                "'_get_AZ' accepts 'str' (particle name) or 'tuple' (A, Z)"
                f", but it received object {type(arg[1])} = {arg[1]}"
            )


class CompositeTarget(object):
    """Definition of composite targets made of multiple (atomic) nuclei.

    Examples of such composite targets are Air, CO_2, HCl, C_2H_60.
    """

    from numpy.random import default_rng

    def __init__(self, component_list, label="", random_state=None):

        if random_state is None:
            self.rng = self.default_rng()
        elif isinstance(random_state, int):
            # Use random_state as a seed
            self.rng = self.default_rng(random_state)
        else:
            # Use random_state as a user provided
            # random number generator
            self.rng = random_state

        self.label = label
        self.ncomponents = 0
        self.component_fractions = []
        self._component_orig_fractions = []
        self.component_A = []
        self.component_Z = []
        self.component_name = []

        for component in component_list:
            if not isinstance(component, tuple):
                raise ValueError("Composite target accepts list of 'tuple's")
            if len(component) == 2:
                self._add_component(component[0], component[1])
            elif len(component) == 3:
                self._add_component(component[0], component[1], component[2])
            else:
                raise ValueError(
                    "'CompositeTarget': wrong component length = "
                    f"{len(component)} for {component}"
                )
        self._normalize()

    def _add_component(self, az, fraction, name=""):
        A, Z = _FromParticleName._get_AZ(az)
        if (name == "") and (isinstance(az, str)):
            self.component_name.append(az)
        else:
            self.component_name.append(name)
        self.component_A.append(A)
        self.component_Z.append(Z)
        self._component_orig_fractions.append(float(fraction))
        self.ncomponents += 1

    def _normalize(self):
        self.component_fractions = np.array(
            [
                f / np.sum(self._component_orig_fractions)
                for f in self._component_orig_fractions
            ]
        )

        self._sort()

        if not (
            len(self.component_name)
            == len(self.component_A)
            == len(self.component_Z)
            == len(self.component_fractions)
            == len(self._component_orig_fractions)
        ):
            raise RuntimeError(
                "CompositeTarget:_normalize len of arrays should be equal"
            )

        if not (np.sum(self.component_fractions) == 1.0):
            raise RuntimeError("Normalization failed")

    def add_component(self, az, fraction, name=""):
        """Add material for composite target.

        Fraction needs relative specification, in percent, number,
        fraction of one, etc. Just make sure it's the same definition
        for all components. Internal list is sorted according to
        relative fractions.
        """
        self._add_component(az, fraction, name)
        self._normalize()

    def _sort(self):
        """Sorts list acording to fraction"""

        def sort_list(k, idcs):
            return [k[i] for i in idcs]

        idcs = np.argsort(self.component_fractions)
        self.component_fractions = self.component_fractions[idcs]
        self._component_orig_fractions = sort_list(self._component_orig_fractions, idcs)
        self.component_A = sort_list(self.component_A, idcs)
        self.component_Z = sort_list(self.component_Z, idcs)
        self.component_name = sort_list(self.component_name, idcs)

    def _get_maximum_AZ(self):
        a_val = self.component_A
        max_ind = a_val.index(max(a_val))
        return self.component_A[max_ind], self.component_Z[max_ind]

    def _get_random_AZ(self):
        """Return randomly an (A, Z) tuple according to the component fraction."""
        ic = self.rng.choice(self.ncomponents, 1, p=self.component_fractions)[0]
        return self.component_A[ic], self.component_Z[ic]

    def __str__(self):
        ostr = "Composite target '" + self.label + "' \n"

        ostr += "Average mass: {0:5.3f}\n".format(
            np.sum([A * f for A, f in zip(self.component_A, self.component_fractions)])
        )

        ostr += "  Nr   |    Name         |  A  |  Z  | fraction\n"
        templ = "  {0:3}  | {1:15s} |{2:4} |{3:4} |  {4:5.4f}\n"
        for i in range(self.ncomponents):
            ostr += templ.format(
                i,
                self.component_name[i],
                self.component_A[i],
                self.component_Z[i],
                self.component_fractions[i],
            )

        return ostr


def _get_particle_input_type(arg):
    if isinstance(arg, int):
        return "pdg_id"
    elif isinstance(arg, str):
        return "string"
    elif _is_AZ_tuple(arg):
        return "tuple"
    elif isinstance(arg, CompositeTarget):
        return "composite_target"
    else:
        raise ValueError(f"Unmaintained parameter type {type(arg)} = {arg}")


def _normalize_particle(particle):
    pdg = None
    nuc_prop = None
    composite_target = None
    received = _get_particle_input_type(particle)

    if received == "pdg_id":
        pdg = particle
    elif received == "tuple":
        nuc_prop = particle
    elif received == "string":
        pdg = _FromParticleName._get_pdg(particle)
    elif received == "composite_target":
        nuc_prop = particle._get_maximum_AZ()
        composite_target = particle

    return pdg, nuc_prop, composite_target


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

    def __init__(
        self,
        *,
        ecm=None,
        plab=None,
        elab=None,
        ekin=None,
        beam=None,
        particle1=None,
        particle2=None,
    ):
        # Catch input errors

        if sum(bool(x) for x in [ecm, plab, elab, ekin, beam]) != 1:
            raise ValueError(
                "Please provide only one of ecm/plab/elab/ekin/beam arguments"
            )

        if particle1:
            p1pdg, nuc1_prop, composite_target = _normalize_particle(particle1)
            if composite_target:
                raise ValueError("Only 2nd parameter could be composite target")
        else:
            raise ValueError(
                '"particle1" parameter is empty. Please provide one of the following'
                " accepted argument types: \nparticle name (str), PDG ID (int), or "
                "nucleus mass & charge (A, Z) tuple"
            )

        if particle2:
            p2pdg, nuc2_prop, self.composite_target = _normalize_particle(particle2)
        else:
            raise ValueError(
                '"particle2" parameter is empty. Please provide one of the following'
                " accepted argument types:\nparticle name (str), PDG ID (int),"
                " nucleus mass & charge (A, Z) tuple, or CompositeTarget object"
            )

        # Store average nucleon mass
        mnuc = 0.5 * (pdata.mass(2212) + pdata.mass(2112))

        # Handle projectile type
        if p1pdg:
            pmass1 = pdata.mass(p1pdg)
            self.A1 = 1
            self.Z1 = pdata.charge(p1pdg)
            self.p1pdg = p1pdg
            self.p1_is_nucleus = False
            info(20, "Particle 1 identified from PDG ID.")
        else:
            pmass1 = mnuc
            self.p1pdg = 2212
            self.A1, self.Z1 = nuc1_prop
            self.p1_is_nucleus = True if self.A1 > 1 else False
            if self.A1 == 1 and self.Z1 == 0:
                self.p1pdg = 2112
            info(20, "Particle 1 is a nucleus.")

        # Handle target type
        if p2pdg:
            pmass2 = pdata.mass(p2pdg)
            self.p2pdg = p2pdg
            self.A2 = 1
            self.Z2 = pdata.charge(p2pdg)
            self.p2_is_nucleus = False
            info(20, "Particle 2 identified from PDG ID.")
        else:
            pmass2 = mnuc
            self.p2pdg = 2212
            self.A2, self.Z2 = nuc2_prop
            self.p2_is_nucleus = True if self.A2 > 1 else False
            if self.A2 == 1 and self.Z2 == 0:
                self.p2pdg = 2112
            info(20, "Particle 2 is a nucleus.")

        info(
            10,
            "Proj. and targ. identified",
            (self.p1pdg, self.A1, self.Z1),
            (self.p2pdg, self.A2, self.Z2),
        )

        # Input specification in center of mass frame
        if ecm:
            self.ecm = ecm
            self.elab = 0.5 * (ecm**2 - pmass1**2 - pmass2**2) / pmass2
            self.plab = self._e2p(self.elab, pmass1)
            info(20, "ecm specified.")
        # Input specification in lab frame
        elif elab:
            if not (elab > pmass1):
                raise RuntimeError("Lab. energy > particle mass required.")
            self.elab = elab
            self.plab = self._e2p(self.elab, pmass1)
            self.ecm = np.sqrt(2.0 * self.elab * pmass2 + pmass2**2 + pmass1**2)
            # self.ecm = np.sqrt((self.elab + pmass2)**2 - self.plab**2)
            info(20, "elab specified.")
        elif ekin:
            self.elab = ekin + pmass1
            self.plab = self._e2p(self.elab, pmass1)
            self.ecm = np.sqrt(2.0 * self.elab * pmass2 + pmass2**2 + pmass1**2)
            # self.ecm = np.sqrt((self.elab + pmass2)**2 - self.plab**2)
            info(20, "ekin specified.")
        elif plab:
            self.plab = plab
            self.elab = np.sqrt(plab**2 + pmass1**2)
            self.ecm = np.sqrt(2.0 * self.elab * pmass2 + pmass2**2 + pmass1**2)
            # self.ecm = np.sqrt((self.elab + pmass2)**2 - self.plab**2)
            info(20, "plab specified.")

        # Input specification as 4-vectors
        elif beam:
            p1, p2 = beam
            s = p1 + p2
            self.ecm = np.sqrt(s[3] ** 2 - np.sum(s[:3] ** 2))
            self.elab = 0.5 * (self.ecm**2 - pmass1**2 + pmass2**2) / pmass2
            self.plab = self._e2p(self.elab, pmass1)
            info(
                20,
                "beam spec: beam, ecms, elab, plab",
                beam,
                self.ecm,
                self.elab,
                self.plab,
            )
        else:
            raise Exception(
                self.__class__.__name__ + "::init(): Define at least ecm or plab"
            )

        self.pmass1 = pmass1
        self.pmass2 = pmass2

        # compute center-of-mass variables
        s = self.ecm**2
        self.s = s
        self.pcm = np.sqrt(
            (s - (pmass1 + pmass2) ** 2) * (s - (pmass1 - pmass2) ** 2)
        ) / (2 * self.ecm)
        self.gamma_cm = (self.elab + pmass2) / self.ecm
        self.betagamma_z_cm = self.plab / self.ecm

        self.boost_def = {
            ("center-of-mass", "laboratory"): self.boost_cms_to_lab,
            ("laboratory", "center-of-mass"): self.boost_lab_to_cms,
        }

        # self.e_range = []

    @property
    def beam_as_4vec(self):
        """Return the projectile target kinematics as 4-vectors. Can be used
        for PHOJET and PYTHIA."""

        p1, p2 = np.array(np.zeros(4), dtype="d"), np.array(np.zeros(4), dtype="d")
        p1[0] = 0.0
        p1[1] = 0.0
        p1[2] = self.pcm
        p1[3] = np.sqrt(self.pcm**2 + self.pmass1**2)
        p2[0] = 0.0
        p2[1] = 0.0
        p2[2] = -self.pcm
        p2[3] = np.sqrt(self.pcm**2 + self.pmass2**2)

        return p1, p2

    def _e2p(self, E, m):
        return np.sqrt((E + m) * (E - m))

    # def __getstate__(self):
    #     if "boost_def" in self.__dict__:
    #         _ = self.__dict__.pop("boost_def")
    #     return self.__dict__

    # def __setstate__(self, state):
    #     self.__dict__ = state
    #     self.boost_def = {
    #         ("center-of-mass", "laboratory"): self.boost_cms_to_lab,
    #         ("laboratory", "center-of-mass"): self.boost_lab_to_cms,
    #     }

    def __ne__(self, other):
        for key, value in other.__dict__.items():
            if key == "boost_def":
                continue
            if value != self.__dict__[key]:
                return True

        return False

    def __eq__(self, other):
        return not self.__ne__(other)

    def __hash__(self):
        return hash(
            "_".join(
                [
                    "{0}_{1}".format(key, self.__dict__[key])
                    for key in sorted(self.__dict__.keys())
                ]
            )
        )

    def __le__(self, other):
        return self.ecm < other.ecm

    def __ge__(self, other):
        return self.ecm > other.ecm

    def __repr__(self):
        ostr = "Event kinematics:\n"
        ostr += "\tecm      : {0:10.5f}\n".format(self.ecm)
        ostr += "\tpcm      : {0:10.5f}\n".format(self.pcm)
        ostr += "\telab     : {0:10.5f}\n".format(self.elab)
        ostr += "\tplab     : {0:10.5f}\n".format(self.plab)
        ostr += "\tgamma_cm : {0:10.5f}\n".format(self.gamma_cm)
        ostr += "\tbgamm_cm : {0:10.5f}\n".format(self.betagamma_z_cm)
        ostr += "\tpdgid 1  : {0:10}\n".format(self.p1pdg)
        ostr += "\tnucprop 1: {0}/{1}\n".format(self.A1, self.Z1)
        ostr += "\tpdgid 2  : {0:10}\n".format(self.p2pdg)
        ostr += "\tnucprop 2: {0}/{1}\n".format(self.A2, self.Z2)

        return ostr

    def apply_boost(self, event, gen_frame, user_frame):
        if (gen_frame, user_frame) not in self.boost_def:
            info(20, "FS boost not applicable", gen_frame, "->", user_frame)
            return
        info(20, "Boosting FS", gen_frame, "->", user_frame)
        self.boost_def[(gen_frame, user_frame)](event)

    def boost_cms_to_lab(self, event):
        """Boosts from center of mass to lab frame.

        Viewed from the target rest frame the center of mass frame
        is moving backwards.
        """
        new_en = self.gamma_cm * event.en + self.betagamma_z_cm * event.pz
        event.pz = self.betagamma_z_cm * event.en + self.gamma_cm * event.pz
        event.en = new_en

    def boost_lab_to_cms(self, event):
        """Boosts from lab to center of mass frame

        Viewed from the target rest frame the center of mass frame
        is moving backwards.
        """
        new_en = self.gamma_cm * event.en - self.betagamma_z_cm * event.pz
        event.pz = -self.betagamma_z_cm * event.en + self.gamma_cm * event.pz
        event.en = new_en


class CenterOfMass(EventKinematics):
    def __init__(self, ecm, particle1, particle2):
        impy_config["user_frame"] = "center-of-mass"
        super().__init__(ecm=ecm, particle1=particle1, particle2=particle2)


class TotalEnergy(float):
    def __new__(cls, value):
        return float.__new__(cls, value)


class KinEnergy(float):
    def __new__(cls, value):
        return float.__new__(cls, value)


class Momentum(float):
    def __new__(cls, value):
        return float.__new__(cls, value)


class FixedTarget(EventKinematics):
    def __init__(self, energy, particle1, particle2):
        impy_config["user_frame"] = "laboratory"

        if isinstance(energy, TotalEnergy):
            super().__init__(elab=energy, particle1=particle1, particle2=particle2)
        elif isinstance(energy, KinEnergy):
            super().__init__(ekin=energy, particle1=particle1, particle2=particle2)
        elif isinstance(energy, Momentum):
            super().__init__(plab=energy, particle1=particle1, particle2=particle2)
        elif isinstance(energy, float):
            super().__init__(elab=energy, particle1=particle1, particle2=particle2)
        else:
            raise ValueError("Please use 'float' (in GeV) when specifying energy")
