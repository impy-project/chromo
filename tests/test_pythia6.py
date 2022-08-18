from impy.kinematics import EventKinematics
from impy.constants import TeV
from impy.models.pythia6 import PYTHIA6Event
from impy.models import Pythia6
from particle import Particle, ParticleNotFound
from numpy.testing import assert_allclose

ekin = EventKinematics(ecm=10 * TeV, p1pdg=2212, p2pdg=2212)
m = Pythia6(ekin, seed=1)
for event in m(1):
    pass


def test_event_type():
    assert isinstance(event, PYTHIA6Event)


def reference_charge(pid):
    try:
        return Particle.from_pdgid(pid).charge
    except ParticleNotFound:
        return 0


def test_charge():
    expected = [reference_charge(pid) for pid in event.id]
    assert_allclose(event.charge, expected)


def test_children():
    assert event.children.shape == (2, len(event))


def test_parents():
    assert event.parents.shape == (2, len(event))
