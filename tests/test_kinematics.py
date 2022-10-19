from impy.kinematics import FixedTarget, TotalEnergy, KinEnergy, Momentum
from impy.constants import GeV
from particle import literals as lp
from pytest import approx


def test_fixed_target():
    ft = FixedTarget(10, "proton", "proton")
    assert ft.elab == 10 * GeV

    ft = FixedTarget(10.0 * GeV, "proton", "proton")
    assert ft.elab == 10 * GeV

    ft = FixedTarget(TotalEnergy(2 * GeV), "proton", "proton")
    assert ft.elab == 2 * GeV

    ft = FixedTarget(KinEnergy(2 * GeV), "proton", "proton")
    et = 2 + (lp.proton.mass / 1e3)
    assert ft.elab == approx(et)

    ft = FixedTarget(Momentum(2 * GeV), "proton", "proton")
    et = (2**2 + (lp.proton.mass / 1e3) ** 2) ** 0.5
    assert ft.elab == approx(et)
