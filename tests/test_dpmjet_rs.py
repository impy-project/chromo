import numpy as np
import pytest

from chromo.kinematics import CenterOfMass, FixedTarget

from .util import baryon_number, charge_number

pytest.importorskip("_dpmjet_rs", reason="DPMJET-rs is not installed")

from chromo.models._extra_models import DpmjetRS


@pytest.fixture(scope="module")
def model():
    return DpmjetRS(CenterOfMass(100, "p", "p"), seed=1)


def _conserved(ev, a, z):
    fs = ev.final_state()
    assert sum(baryon_number(p) for p in fs.pid) == a
    assert sum(charge_number(p) for p in fs.pid) == z
    return fs


def test_proton_proton(model):
    model.kinematics = CenterOfMass(100, "p", "p")
    for ev in model(50):
        fs = _conserved(ev, 2, 2)
        assert fs.en.sum() == pytest.approx(100, rel=1e-6)
        assert np.abs(fs.pz.sum()) < 1e-4
    cs = model.cross_section()
    assert cs.total == pytest.approx(cs.inelastic + cs.elastic)
    assert 35 < cs.inelastic < 45
    assert cs.diffractive_sum < cs.inelastic


def test_nuclei(model):
    # residual nuclei carry the spectators: baryon number and charge are
    # conserved; the energy up to their binding and Fermi motion
    model.kinematics = CenterOfMass(50, (12, 6), (14, 7))
    for ev in model(20):
        fs = _conserved(ev, 26, 13)
        assert fs.en.sum() == pytest.approx(26 * 25, rel=0.02)
    assert 800 < model.cross_section().prod < 1100


def test_pion_fixed_target(model):
    model.kinematics = FixedTarget(100, "pi-", (12, 6))
    for ev in model(20):
        _conserved(ev, 12, 5)
    assert 150 < model.cross_section().prod < 250


def test_reproducible():
    kin = CenterOfMass(200, "p", (14, 7))
    a = DpmjetRS(kin, seed=7)
    b = DpmjetRS(kin, seed=7)
    for ea, eb in zip(a(5), b(5)):
        assert np.array_equal(ea.pid, eb.pid)
        assert np.array_equal(ea.px, eb.px)


def test_unsupported(model):
    with pytest.raises(ValueError):
        model.kinematics = CenterOfMass(100, "gamma", "p")


def test_decays_before_cascade(model):
    # resonances decay inside the model; set_stable keeps a species
    model.kinematics = CenterOfMass(200, "p", (12, 6))
    pids = np.concatenate([ev.final_state().pid for ev in model(50)])
    assert not np.isin(np.abs(pids), [113, 213, 223, 2224, 3114]).any()
    model.set_stable(223, True)
    pids = np.concatenate([ev.final_state().pid for ev in model(50)])
    assert np.isin(pids, [223]).any()
    model.set_stable(223, False)
