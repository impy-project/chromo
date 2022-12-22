from impy.kinematics import (
    FixedTarget,
    TotalEnergy,
    KinEnergy,
    Momentum,
    EventFrame,
    CompositeTarget,
)
from impy.constants import GeV, MeV
from particle import literals as lp
from pytest import approx


def test_fixed_target():
    ft = FixedTarget(10, "proton", "proton")
    assert ft.elab == 10 * GeV
    assert ft.frame == EventFrame.FIXED_TARGET

    ft = FixedTarget(10.0 * GeV, "proton", "proton")
    assert ft.elab == 10 * GeV
    assert ft.frame == EventFrame.FIXED_TARGET

    ft = FixedTarget(TotalEnergy(2 * GeV), "proton", "proton")
    assert ft.elab == 2 * GeV
    assert ft.frame == EventFrame.FIXED_TARGET

    ft = FixedTarget(KinEnergy(2 * GeV), "proton", "proton")
    et = 2 + (lp.proton.mass * MeV)
    assert ft.elab == approx(et)
    assert ft.frame == EventFrame.FIXED_TARGET

    ft = FixedTarget(Momentum(2 * GeV), "proton", "proton")
    et = (2**2 + (lp.proton.mass * MeV) ** 2) ** 0.5
    assert ft.elab == approx(et)
    assert ft.frame == EventFrame.FIXED_TARGET


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
