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

# I'm used to global configs, but I understand that this is often not good.
# It's a must, though, to have one global instance of the particle tables
# accessing a locally cached database. The XML parsing time was often a
# problem in the past and this is a serious performance concern.

# from somewhere_in_impy import pdata
# For now:
import numpy as np
from impy.common import pdata, impy_config
from impy.util import info

class CompositeTarget(object):
    """Definition of composite targets made of multiple (atomic) nuclei.

    Examples of such composite targets are Air, CO_2, HCl, C_2H_60.
    """

    def __init__(self, label=''):
        self.label = label
        self.ncomponents = 0
        self.component_fractions = []
        self._component_orig_fractions = []
        self.component_A = []
        self.component_Z = []
        self.component_name = []

    def add_component(self, name, A, Z, fraction):
        """Add material for composite target.
        
        Fraction needs relative specification, in percent, number,
        fraction of one, etc. Just make sure it's the same definition
        for all components. Internal list is sorted according to
        relative fractions.
        """

        self.component_name.append(name)
        self.component_A.append(A)
        self.component_Z.append(Z)
        self._component_orig_fractions.append(float(fraction))
        self.component_fractions = np.array([
            f / np.sum(self._component_orig_fractions)
            for f in self._component_orig_fractions
        ])
        self.ncomponents += 1

        self._sort()

        assert (len(self.component_name) == len(self.component_A) == len(
            self.component_Z) == len(self.component_fractions) == len(
                self._component_orig_fractions))
        assert np.sum(self.component_fractions) == 1.

    def _sort(self):
        """Sorts list acording to fraction"""
        sort_list = lambda l, idcs: [l[i] for i in idcs]
        idcs = np.argsort(self.component_fractions)
        self.component_fractions = self.component_fractions[idcs]
        self._component_orig_fractions = sort_list(
            self._component_orig_fractions, idcs)
        self.component_A = sort_list(self.component_A, idcs)
        self.component_Z = sort_list(self.component_Z, idcs)
        self.component_name = sort_list(self.component_name, idcs)

    def get_random_AZ(self):
        """Return randomly an (A, Z) tuple according to the component fraction."""
        from numpy.random import choice

        ic = choice(self.ncomponents, 1, p=self.component_fractions)
        return self.component_A[ic], self.component_Z[ic]

    def __repr__(self):
        ostr = "Composite target '" + self.label + "' \n"

        ostr += "Average mass: {0:4.1f}\n".format(
            np.sum([
                A * f
                for A, f in zip(self.component_A, self.component_fractions)
            ]))

        ostr += "  Nr   |    Name    |  A  |  Z  | fraction\n"
        templ = "  {0:3}  | {1:10s} |{2:4} |{3:4} |  {4:3.2f}\n"
        for i in range(self.ncomponents):
            ostr += templ.format(i, self.component_name[i],
                                 self.component_A[i], self.component_Z[i],
                                 self.component_fractions[i])

        return ostr


class EventKinematics(object):
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
        (float) p1pdg    : PDG ID of the projectile
        (float) p2pdg    : PDG ID of the target
        (tuple) beam     : Specification as tuple of 4-vectors (np.array)s
        (tuple) nuc1prop : Projectile nucleus mass & charge (A, Z)
        (tuple) nuc2prop : Target nucleus mass & charge (A, Z)
     
    """

    def __init__(self,
                 ecm=None,
                 plab=None,
                 p1pdg=None,
                 p2pdg=None,
                 beam=None,
                 nuc1_prop=None,
                 nuc2_prop=None):

        # Catch input errors
        if beam is not None:
            raise NotImplementedError('Not properly implmented or not needed')

        if ecm and plab:
            raise Exception(self.__class__.__name__ +
                            '::init(): Define either ecm or plab')

        if p1pdg and nuc1_prop:
            raise Exception(self.__class__.__name__ +
                            '::init(): Define either particle id or ' +
                            'nuclear protperties for side 1.')

        if p2pdg and nuc2_prop:
            raise Exception(self.__class__.__name__ +
                            '::init(): Define either particle id or ' +
                            'nuclear protperties for side 2.')

        # Store average nucleon mass
        mnuc = 0.5 * (pdata.mass(2212) + pdata.mass(2112))

        # Handle projectile type
        if p1pdg:
            pmass1 = pdata.mass(p1pdg)
            self.A1, self.Z1 = 1, 1
            self.p1pdg = p1pdg
        else:
            pmass1 = mnuc
            self.p1pdg = 2212
            self.A1, self.Z1 = nuc1_prop

        # Handle target type
        if p2pdg:
            try:
                pmass2 = pdata.mass(p2pdg)
            except KeyError:
                pmass2 = pdata.mass(p1pdg)

            self.p2pdg = p2pdg
            if p2pdg > 0:
                self.A2, self.Z2 = 1, 1
            else:
                self.A2, self.Z2 = 1, -1
        else:
            pmass2 = mnuc
            self.p2pdg = 2212
            self.A2, self.Z2 = nuc2_prop
        info(10, 'Proj. and targ. identified', (self.p1pdg, self.A1, self.Z1),
             (self.p2pdg, self.A2, self.Z2))

        # Input specification in center of mass frame
        if ecm and not (plab or beam):
            self.ecm = ecm
            self.elab = 0.5 * (ecm**2 - pmass1**2 + pmass2**2) / pmass2
            self.plab = self._e2p(self.elab, pmass1)
            info(20, 'Center of mass variables specified.')

        # Input specification in lab frame
        elif plab and not (ecm or beam):
            self.plab = plab
            self.elab = np.sqrt(plab**2 + pmass1**2)
            self.ecm = np.sqrt((self.elab + pmass2)**2 - self.plab**2)
            info(20, 'Lab variables specified.')

        # Input specification as 4-vectors
        elif beam and not (ecm or plab):
            # self.beam = beam
            # p1 = np.array([0, 0, self._e2p(beam[0][2], pmass1), beam[0][3]])
            # p2 = np.array([0, 0, self._e2p(beam[1][2], pmass2), beam[0][3]])
            p1, p2 = beam
            s = p1 + p2
            self.ecm = np.sqrt(s[3]**2 - np.sum(s[:3]**2))
            self.elab = 0.5 * (self.ecm**2 - pmass1**2 + pmass2**2) / pmass2
            self.plab = self._e2p(self.elab, pmass1)
            info(20, 'beam spec: beam, ecms, elab, plab', beam, self.ecm,
                 self.elab, self.plab)
        else:
            raise Exception(self.__class__.__name__ +
                            '::init(): Define at least ecm or plab')

        self.pmass1 = pmass1
        self.pmass2 = pmass2

        # compute center-of-mass variables
        s = self.ecm**2
        self.s = s
        self.pcm = np.sqrt(s**2 - 2 * (s * (pmass1 + pmass2) + pmass1 * pmass2)
                           + pmass1**2 + pmass2**2) / (2 * self.ecm)
        self.gamma_cm = (self.elab + pmass2) / self.ecm
        self.betagamma_z_cm = self.plab / self.ecm

        self.boost_def = {
            ("center-of-mass", "laboratory") : self.boost_cms_to_lab,
            ("laboratory", "center-of-mass") : self.boost_lab_to_cms
        }

        # self.e_range = []

    def _e2p(self, E, m):
        return np.sqrt((E + m) * (E - m))

    def __ne__(self, other):
        return self.__dict__ != other.__dict__

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __hash__(self):
        return hash('_'.join([
            '{0}_{1}'.format(key, self.__dict__[key])
            for key in sorted(self.__dict__.keys())
        ]))

    def __le__(self, other):
        return self.ecm < other.ecm

    def __ge__(self, other):
        return self.ecm > other.ecm

    def __repr__(self):
        ostr = 'Event kinematics:\n'
        ostr += '\tecm      : {0:10.5f}\n'.format(self.ecm)
        ostr += '\tpcm      : {0:10.5f}\n'.format(self.pcm)
        ostr += '\telab     : {0:10.5f}\n'.format(self.elab)
        ostr += '\tplab     : {0:10.5f}\n'.format(self.plab)
        ostr += '\tgamma_cm : {0:10.5f}\n'.format(self.gamma_cm)
        ostr += '\tbgamm_cm : {0:10.5f}\n'.format(self.betagamma_z_cm)
        ostr += '\tpdgid 1  : {0:10}\n'.format(self.p1pdg)
        ostr += '\tnucprop 1: {0}/{1}\n'.format(self.A1, self.Z1)
        ostr += '\tpdgid 2  : {0:10}\n'.format(self.p2pdg)
        ostr += '\tnucprop 2: {0}/{1}\n'.format(self.A2, self.Z2)

        return ostr

    def apply_boost(self, event, gen_frame, user_frame):
        if (gen_frame, user_frame) not in self.boost_def:
            info(20, 'FS boost not applicable',gen_frame, '->', user_frame)
            return
        info(20, 'Boosting FS', gen_frame, '->', user_frame)
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