import numpy as np
import pytest
from particle import literals as lp
from pytest import approx

from chromo.constants import nucleon_mass
from chromo.kinematics import (
    CenterOfMass,
    CompositeTarget,
    EventFrame,
    EventKinematicsMassless,
    EventKinematicsWithRestframe,
    FixedTarget,
    GeV,
    KinEnergy,
    MeV,
    Momentum,
    TotalEnergy,
)
from chromo.util import (
    AZ2pdg,
    elab2ecm,
    energy2momentum,
    mass,
    momentum2energy,
)


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

    with pytest.raises(TypeError):
        FixedTarget(100 * GeV, t, "p")


def test_copy():
    target = CompositeTarget([("N", 0.78), ("O", 0.21), ("Ar", 0.01)], label="air")
    a = CenterOfMass(10, "p", target)
    b = a.copy()
    assert a == b


def test_kinematics_init_ecm():
    # Test initialization with ecm argument
    k = EventKinematicsWithRestframe("proton", "neutron", ecm=10)
    assert k.frame == EventFrame.CENTER_OF_MASS
    assert k.ecm == 10
    assert k.elab == approx(52.2778, rel=1e-3)
    assert k.ekin == approx(k.elab - k.m1, rel=1e-3)
    assert k.plab == approx(energy2momentum(k.elab, k.m1), rel=1e-3)


def test_kinematics_init_beam():
    # Test initialization with beam argument
    k = EventKinematicsWithRestframe(
        "proton", "neutron", beam=(10.0, -4.0), frame=EventFrame.CENTER_OF_MASS
    )
    assert k.ecm == approx(12.818, rel=1e-3)
    k_ref = EventKinematicsWithRestframe("proton", "neutron", ecm=k.ecm)
    assert k == k_ref


def test_kinematics_init_elab():
    # Test initialization with elab argument
    k = EventKinematicsWithRestframe("proton", "neutron", elab=15.0)
    assert k.frame == EventFrame.FIXED_TARGET
    k_ref = EventKinematicsWithRestframe(
        "proton", "neutron", ecm=k.ecm, frame=EventFrame.FIXED_TARGET
    )
    assert k == k_ref


def test_kinematics_init_ekin():
    # Test initialization with ekin argument
    k = EventKinematicsWithRestframe("proton", "neutron", ekin=8)
    assert k.frame == EventFrame.FIXED_TARGET
    assert k.ecm == approx(
        elab2ecm(8 + nucleon_mass, nucleon_mass, nucleon_mass), rel=1e-3
    )
    assert k.elab == approx(8 + nucleon_mass, rel=1e-3)
    assert k.ekin == 8
    assert k.plab == approx(energy2momentum(8 + nucleon_mass, nucleon_mass), rel=1e-3)


def test_kinematics_init_invalid_input():
    # Test initialization with invalid input
    with pytest.raises(ValueError):
        EventKinematicsWithRestframe("proton", None, ecm=10, plab=5)

    with pytest.raises(ValueError):
        EventKinematicsWithRestframe(None, "neutron", ecm=10, plab=5)

    with pytest.raises(ValueError):
        EventKinematicsWithRestframe("proton", "neutron", ecm=10, plab=5, elab=15)

    with pytest.raises(ValueError):
        EventKinematicsWithRestframe("proton", "neutron", ecm=10, plab=5, ekin=8)

    with pytest.raises(ValueError):
        EventKinematicsWithRestframe("proton", "neutron", ecm=10, plab=5, beam=(5, 3))

    with pytest.raises(ValueError):
        EventKinematicsWithRestframe("photon", "photon", ecm=10.0)

    with pytest.raises(ValueError):
        EventKinematicsWithRestframe("photon", "e+", virtuality=(0.7, 0.8))


def test_kinematics_virtuality():
    # Test initialization with virtuality argument
    k = EventKinematicsMassless("photon", "photon", ecm=10.0, virtuality=(0.5, 0.3))
    assert k.virt_p1 == 0.5
    assert k.virt_p2 == 0.3

    k = EventKinematicsWithRestframe("photon", "e+", ecm=10.0, virtuality=0.7)
    assert k.virt_p1 == 0.7
    assert k.virt_p2 == 0.0


def test_kinematics_composite_target():
    # Test initialization with CompositeTarget
    target = CompositeTarget([("N", 3), ("O", 1)])
    with pytest.raises(TypeError):
        EventKinematicsWithRestframe(target, "neutron", ecm=10)

    k = EventKinematicsWithRestframe(
        "proton", target, ecm=10, frame=EventFrame.CENTER_OF_MASS
    )
    assert k.frame == EventFrame.CENTER_OF_MASS
    assert k.ecm == 10
    k_ref = EventKinematicsWithRestframe("proton", target, ecm=10)
    assert k == k_ref


def test_kinematics_beam_data():
    # Test beam data
    k = EventKinematicsWithRestframe("proton", "neutron", beam=(5, 3))
    assert k.beams[0][2] == 5
    assert k.beams[1][2] == 3
    assert k.beams[0][3] == approx(momentum2energy(5, nucleon_mass), rel=1e-3)
    assert k.beams[1][3] == approx(momentum2energy(3, nucleon_mass), rel=1e-3)


def test_kinematics_gamma_cm():
    # Test gamma_cm calculation
    k = EventKinematicsWithRestframe("proton", "neutron", ecm=10)
    assert k._gamma_cm == (k.elab + k.m2) / k.ecm


def test_kinematics_betagamma_cm():
    # Test betagamma_cm calculation
    k = EventKinematicsWithRestframe("proton", "neutron", ecm=10)
    assert k._betagamma_cm == approx(k.plab / k.ecm, rel=1e-3)


def test_kinematics_m1_m2():
    # Test m1 and m2 values
    k = EventKinematicsWithRestframe("proton", "neutron", ecm=10)
    assert np.allclose(k.m1, mass(2212))
    assert np.allclose(k.m2, mass(2112))
