from chromo.kinematics import (
    CenterOfMass,
    FixedTarget,
    TotalEnergy,
    KinEnergy,
    Momentum,
    EventFrame,
    CompositeTarget,
    GeV,
    MeV,
)
from chromo.constants import nucleon_mass
from chromo.util import AZ2pdg, energy2momentum
from particle import literals as lp
from pytest import approx
import pytest
import numpy as np


def test_CompositeTarget_repr():
    t = CompositeTarget([("N", 3), ("O", 1)])
    assert t.A == 16
    assert t.Z == 8
    assert t.components == (1000070140, 1000080160)
    assert int(t) == int(t.components[1])
    assert abs(t) == int(t.components[1])
    assert repr(t) == "CompositeTarget([('N14', 0.75), ('O16', 0.25)])"

    t = CompositeTarget([("N", 3), ("O", 1)], label="air")
    assert repr(t) == "CompositeTarget([('N14', 0.75), ('O16', 0.25)], label='air')"


def test_CompositeTarget_copy():
    target = CompositeTarget([("N", 0.78), ("O", 0.21), ("Ar", 0.01)], label="air")
    target1 = target.copy()
    assert target1 is not target
    assert target1 == target
    assert target != "p"


def test_fixed_target():
    x = 2 * GeV

    ft = FixedTarget(TotalEnergy(x), "proton", "proton")
    assert ft.plab < x
    assert ft.elab == x
    assert ft.frame == EventFrame.FIXED_TARGET
    # default is to interpret x as total energy
    assert ft == FixedTarget(x, "proton", "proton")

    ft = FixedTarget(KinEnergy(x), "proton", "proton")
    et = x + (lp.proton.mass * MeV)
    assert ft.elab == approx(et, rel=1e-3)
    assert ft.frame == EventFrame.FIXED_TARGET

    ft = FixedTarget(Momentum(x), "proton", "proton")
    et = (x**2 + (lp.proton.mass * MeV) ** 2) ** 0.5
    assert ft.plab == x
    assert ft.elab > x
    assert ft.elab == approx(et, rel=1e-3)
    assert ft.frame == EventFrame.FIXED_TARGET

    ft = FixedTarget(x, "proton", "He")
    assert ft.p1 == lp.proton.pdgid
    assert ft.p2 == AZ2pdg(4, 2)
    # check that ecm is in nucleon-nucleon collision system
    p1 = np.array([energy2momentum(x, lp.proton.mass * MeV), x])
    p2 = np.array([0, nucleon_mass])
    ps = p1 + p2
    ecm = (ps[1] ** 2 - ps[0] ** 2) ** 0.5
    assert ft.ecm == approx(ecm, rel=1e-3)

    x = 32 * GeV
    ft = FixedTarget(x, "O", "He")
    assert ft.p1 == AZ2pdg(16, 8)
    assert ft.p2 == AZ2pdg(4, 2)
    # check that ecm is in nucleon-nucleon collision system
    p1 = np.array([energy2momentum(x, nucleon_mass), x])
    p2 = np.array([0, nucleon_mass])
    ps = p1 + p2
    ecm = (ps[1] ** 2 - ps[0] ** 2) ** 0.5
    assert ft.ecm == approx(ecm, rel=1e-3)


def test_fixed_target_bad_input():
    with pytest.raises(ValueError):
        FixedTarget(0.1 * GeV, "p", "p")

    t = CompositeTarget([("N", 3), ("O", 1)])

    with pytest.raises(ValueError):
        FixedTarget(100 * GeV, t, "p")


def test_copy():
    target = CompositeTarget([("N", 0.78), ("O", 0.21), ("Ar", 0.01)], label="air")
    a = CenterOfMass(10, "p", target)
    b = a.copy()
    assert a == b
