from impy.models import Pythia6
from impy.kinematics import EventKinematics
from impy.constants import TeV
from particle import Particle
from numpy.testing import assert_allclose

ekin = EventKinematics(ecm=10 * TeV, p1pdg=2212, p2pdg=2212)
pythia6 = Pythia6(ekin)


def test_charge():
    for event in pythia6(1):
        pass

    expected = [Particle.from_pdgid(pid).charge for pid in event.p_ids]
    assert_allclose(event.charge, expected)
