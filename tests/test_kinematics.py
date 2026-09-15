from types import SimpleNamespace

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
    boost_event,
    boost_vector,
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


def test_boost_event_analytic_z():
    # particle at rest, E=2, boost with b=(0,0,0.6): gamma=1.25
    ev = SimpleNamespace(
        en=np.array([2.0]), px=np.array([0.0]), py=np.array([0.0]), pz=np.array([0.0])
    )
    boost_event(ev, (0, 0, 0.6))
    assert ev.en[0] == approx(2.5)
    assert ev.pz[0] == approx(-1.5)
    assert ev.px[0] == 0
    assert ev.py[0] == 0
    # b=(0,0,1) should be rejected
    with pytest.raises(ValueError):
        boost_event(ev, (0, 0, 1.0))
    # zero boost is a no-op
    ev2 = SimpleNamespace(
        en=np.array([1.0]), px=np.array([0.1]), py=np.array([0.2]), pz=np.array([0.3])
    )
    boost_event(ev2, (0, 0, 0))
    assert ev2.en[0] == 1.0


def test_boost_event_generic_preserves_invariants():
    rng = np.random.default_rng(1)
    m = rng.uniform(0.1, 5, 100)
    px, py, pz = (rng.normal(0, 3, 100) for _ in range(3))
    en = np.sqrt(m**2 + px**2 + py**2 + pz**2)
    ev = SimpleNamespace(en=en.copy(), px=px.copy(), py=py.copy(), pz=pz.copy())
    b = (0.1, -0.2, 0.35)
    boost_event(ev, b)
    inv2 = ev.en**2 - ev.px**2 - ev.py**2 - ev.pz**2
    assert inv2 == approx(m**2, rel=1e-10)
    # total four-momentum transforms like a single four-vector
    p_from = np.array([px.sum(), py.sum(), pz.sum(), en.sum()])
    p_to = np.array([ev.px.sum(), ev.py.sum(), ev.pz.sum(), ev.en.sum()])
    assert boost_vector(p_from, p_to) == approx(np.array(b))
    # inverse boost restores the original four-vectors
    boost_event(ev, -np.array(b))
    assert ev.en == approx(en)
    assert ev.px == approx(px)
    assert ev.py == approx(py)
    assert ev.pz == approx(pz)


def test_boost_vector_roundtrip():
    rng = np.random.default_rng(7)
    for _ in range(100):
        m = rng.uniform(0.1, 5)
        p = rng.normal(0, 3, 3)
        P = np.array([*p, np.sqrt(m**2 + p @ p)])
        b = rng.uniform(-0.9, 0.9, 3)
        while b @ b > 0.98:
            b = rng.uniform(-0.9, 0.9, 3)
        ev = SimpleNamespace(
            en=np.array([P[3]]),
            px=np.array([P[0]]),
            py=np.array([P[1]]),
            pz=np.array([P[2]]),
        )
        boost_event(ev, b)
        Pp = np.array([ev.px[0], ev.py[0], ev.pz[0], ev.en[0]])
        assert boost_vector(P, Pp) == approx(b, abs=1e-10)


def test_apply_boost_cms2ft_matches_old_collinear():
    k = EventKinematicsWithRestframe("proton", "neutron", elab=1000)
    ev = SimpleNamespace(
        en=np.array([1.0, 3.0]),
        px=np.array([0.1, -2.0]),
        py=np.array([0.5, 1.0]),
        pz=np.array([0.2, 4.0]),
    )
    en0, pz0 = ev.en.copy(), ev.pz.copy()
    k.apply_boost(ev, EventFrame.CENTER_OF_MASS)
    g, bg = k._gamma_cm, k._betagamma_cm
    assert ev.en == approx(g * en0 + bg * pz0, rel=1e-12)
    assert ev.pz == approx(bg * en0 + g * pz0, rel=1e-12)
    k.apply_boost(ev, EventFrame.CENTER_OF_MASS, inverse=True)
    assert ev.en == approx(en0)
    assert ev.pz == approx(pz0)
    # UHECR energies: boost must stay exact in (gamma, betagamma),
    # a boost reconstructed from b alone loses ~gamma**2 * eps precision
    k_uhe = EventKinematicsWithRestframe("proton", "neutron", elab=1e11)
    ev_uhe = SimpleNamespace(
        en=np.array([1.0]), px=np.array([0.1]), py=np.array([0.5]), pz=np.array([0.2])
    )
    k_uhe.apply_boost(ev_uhe, EventFrame.CENTER_OF_MASS)
    g, bg = k_uhe._gamma_cm, k_uhe._betagamma_cm
    assert ev_uhe.en[0] == approx(g + bg * 0.2, rel=1e-10)
    assert ev_uhe.pz[0] == approx(bg + g * 0.2, rel=1e-10)


def test_apply_boost_generic_frame_pA():
    # asymmetric collider configuration from issue #182
    e_beam = 6.8e3
    k = EventKinematicsWithRestframe("p", "O", beam=(e_beam, -e_beam * 8 / 16))
    assert k.frame == EventFrame.GENERIC
    total_generic = k.beams[0] + k.beams[1]
    # a system at rest in the CMS must be boosted to the total beam momentum
    ev = SimpleNamespace(
        en=np.array([k.ecm / 2, k.ecm / 2]),
        px=np.array([0.0, 0.0]),
        py=np.array([0.0, 0.0]),
        pz=np.array([0.0, 0.0]),
    )
    k.apply_boost(ev, EventFrame.CENTER_OF_MASS)
    total = np.array([ev.px.sum(), ev.py.sum(), ev.pz.sum(), ev.en.sum()])
    assert total == approx(total_generic, rel=1e-10)
    # transverse momenta unchanged by a boost along the beam axis
    assert ev.px.sum() == approx(0, abs=1e-12)
    assert ev.py.sum() == approx(0, abs=1e-12)
    # a generic frame specified as fixed target matches FixedTarget
    kft = EventKinematicsWithRestframe("proton", "neutron", elab=1000)
    kgen = EventKinematicsWithRestframe("proton", "neutron", beam=(kft.plab, 0))
    assert kgen.frame == EventFrame.GENERIC
    ev1 = SimpleNamespace(
        en=np.array([1.0, 3.0]),
        px=np.array([0.1, -2.0]),
        py=np.array([0.5, 1.0]),
        pz=np.array([0.2, 4.0]),
    )
    ev2 = SimpleNamespace(
        en=ev1.en.copy(), px=ev1.px.copy(), py=ev1.py.copy(), pz=ev1.pz.copy()
    )
    kft.apply_boost(ev1, EventFrame.CENTER_OF_MASS)
    kgen.apply_boost(ev2, EventFrame.CENTER_OF_MASS)
    assert ev2.en == approx(ev1.en, rel=1e-10)
    assert ev2.pz == approx(ev1.pz, rel=1e-10)


def test_apply_boost_generic_symmetric_is_cms():
    k = EventKinematicsWithRestframe("proton", "neutron", beam=(5, -5))
    assert k.frame == EventFrame.GENERIC
    ev = SimpleNamespace(
        en=np.array([1.0]),
        px=np.array([0.1]),
        py=np.array([0.5]),
        pz=np.array([0.2]),
    )
    k.apply_boost(ev, EventFrame.CENTER_OF_MASS)
    assert ev.en[0] == approx(1.0, abs=1e-12)
    assert ev.px[0] == approx(0.1, abs=1e-12)
    assert ev.py[0] == approx(0.5, abs=1e-12)
    assert ev.pz[0] == approx(0.2, abs=1e-12)
    # boosting from the generic frame is not supported
    kft = EventKinematicsWithRestframe("proton", "neutron", elab=1000)
    with pytest.raises(NotImplementedError):
        kft.apply_boost(ev, EventFrame.GENERIC)


def test_apply_boost_generic_inverse_roundtrip():
    k = EventKinematicsWithRestframe("p", "O", beam=(6.8e3, -3.4e3))
    ev = SimpleNamespace(
        en=np.array([1.0, 3.0]),
        px=np.array([0.1, -2.0]),
        py=np.array([0.5, 1.0]),
        pz=np.array([0.2, 4.0]),
    )
    ref = (ev.en.copy(), ev.px.copy(), ev.py.copy(), ev.pz.copy())
    k.apply_boost(ev, EventFrame.CENTER_OF_MASS)
    assert ev.en**2 - ev.px**2 - ev.py**2 - ev.pz**2 == approx(
        ref[0] ** 2 - ref[1] ** 2 - ref[2] ** 2 - ref[3] ** 2, rel=1e-10
    )
    k.apply_boost(ev, EventFrame.CENTER_OF_MASS, inverse=True)
    assert ev.en == approx(ref[0])
    assert ev.px == approx(ref[1])
    assert ev.py == approx(ref[2])
    assert ev.pz == approx(ref[3])
