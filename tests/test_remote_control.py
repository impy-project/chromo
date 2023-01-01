from impy.models import Phojet193
from impy.kinematics import CenterOfMass


def test_phojet_1():
    model = Phojet193(seed=1)

    kin = CenterOfMass(100, "p", "p")

    events1 = []
    for event in model(kin, 10):
        events1.append(event)

    assert len(events1) == 10

    for event in events1:
        assert len(event) > 0

    events2 = []
    for event in model(kin, 5):
        events2.append(event)
    assert len(events2) == 5

    assert events2 != events1[:5]


def test_phojet_2():
    model = Phojet193(seed=1)

    kin = CenterOfMass(100, "p", "p")

    events1 = []
    for event in model(kin, 10):
        events1.append(event)

    assert len(events1) == 10

    for event in events1:
        assert len(event) > 0

    kin = CenterOfMass(100, "pi-", "p")

    events2 = []
    for event in model(kin, 5):
        events2.append(event)
    assert len(events2) == 5

    assert events2 != events1[:5]


def test_phojet_3():
    model = Phojet193(seed=1)

    prev = 0
    for en in (10, 100, 1000):
        kin = CenterOfMass(en, "p", "p")
        c = model.cross_section(kin)
        assert c.inelastic > 10
        assert c.inelastic > prev
        prev = c.inelastic


def test_phojet_4():
    model = Phojet193(seed=1)

    # this calls an attribute of the original class
    assert len(model.projectiles) > 0
