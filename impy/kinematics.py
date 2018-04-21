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

class EventKinematics():
    def __init__(self,
                 ecm=None,
                 plab=None,
                 p1pdg=None,
                 p2pdg=None,
                 beam=None,
                 nuc1_prop=None,
                 nuc2_prop=None):

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

        masses = {
            2212: 0.93827,
            2112: 0.93957,
            211: 0.13957,
            321: 0.49368,
            3122: 1.11568,
            22: 0.
        }
        # In case the list is not sufficient load ParticleDataTool
        # on demand
        self.pdata = None
        # Store masses of particles on side 1 and 2
        pmass1, pmass2 = None, None
        # Store average nucleon mass
        mnuc = 0.5 * (masses[2212] + masses[2112])

        if p1pdg:
            try:
                pmass1 = masses[p1pdg]
            except KeyError:
                pmass1 = pdata.mass(p1pdg)
            self.A1, self.Z1 = 1, 1
            self.p1pdg = p1pdg
        else:
            pmass1 = mnuc
            self.p1pdg = 2212
            self.A1, self.Z1 = nuc1_prop

        if p2pdg:
            try:
                pmass2 = masses[p2pdg]
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

        if ecm and not (plab or beam):
            self.ecm = ecm
            self.elab = 0.5 * (ecm**2 - pmass1**2 + pmass2**2) / pmass2
            self.plab = self.e2p(self.elab, pmass1)

        elif plab and not (ecm or beam):
            self.plab = plab
            self.elab = np.sqrt(plab**2 + pmass1**2)
            self.ecm = np.sqrt((self.elab + pmass2)**2 - self.plab**2)

        elif beam and not (ecm or plab):
            self.beam = beam
            p1 = np.array([0, 0, self.e2p(beam[0], pmass1), beam[0]])
            p2 = np.array([0, 0, -self.e2p(beam[1], pmass1), beam[1]])
            s = p1 + p2
            self.ecm = np.sqrt(s[3]**2 - np.sum(s[:3]**2))
            self.elab = 0.5 * (self.ecm**2 - pmass1**2 + pmass2**2) / pmass2
            self.plab = self.e2p(self.elab, pmass1)
            print('EventKinematics: beam, ecms, elab, plab', self.beam,
                  self.ecm, self.elab, self.plab)
        else:
            raise Exception(self.__class__.__name__ +
                            '::init(): Define at least ecm or plab')

        self.pmass1 = pmass1
        self.pmass2 = pmass2

        # compute pcm
        s = self.ecm**2
        self.s = s
        self.pcm = np.sqrt(s**2 - 2 * (s * (pmass1 + pmass2) + pmass1 * pmass2)
                           + pmass1**2 + pmass2**2) / (2 * self.ecm)
        self.gamma_cm = (self.elab + pmass2) / self.ecm
        self.betagamma_cm = self.plab / self.ecm

        self.e_range = []

    def e2p(self, E, m):
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
        ostr += '\tbgamm_cm : {0:10.5f}\n'.format(self.betagamma_cm)
        ostr += '\tpdgid 1  : {0:10}\n'.format(self.p1pdg)
        ostr += '\tnucprop 1: {0}/{1}\n'.format(self.A1, self.Z1)
        ostr += '\tpdgid 2  : {0:10}\n'.format(self.p2pdg)
        ostr += '\tnucprop 2: {0}/{1}\n'.format(self.A2, self.Z2)

        return ostr
        