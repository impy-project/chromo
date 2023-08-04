import chromo
import numpy as np
import pytest
import sys

from chromo.decay_handler import Pythia8DecayHandler
from chromo.util import get_all_models
from .util import run_in_separate_process

if sys.platform == "win32":
    pytest.skip("Pythia8 does not run on windows", allow_module_level=True)


def run_decay_handler(model, evt_kin, stable_particles, decaying_particles):
    seed = 1
    decay_handler = Pythia8DecayHandler(
        seed=seed,
        extra_stable_pids=stable_particles,
        extra_decaying_pids=decaying_particles,
    )
    gen = model(evt_kin, seed=seed)

    for event in gen(10):
        event0 = event.copy()
        decay_handler(event)

        # Assert that decay has happend
        event_dec = event[0 : len(event0)]
        assert np.sum(event_dec.status != 1) >= np.sum(event0.status != 1)

        # Assert that all particles that should decay have been decayed
        not_decayed0 = np.sum(
            np.isin(event0.pid, decay_handler.all_decaying_pids) & (event0.status == 1)
        )
        not_decayed1 = np.sum(
            np.isin(event.pid, decay_handler.all_decaying_pids) & (event.status == 1)
        )
        if not_decayed0 > 0:
            assert not_decayed1 == 0

        # Assert that stable particles haven't decayed
        stable0 = np.isin(event0.pid, decay_handler.all_stable_pids) & (
            event0.status == 1
        )
        stable1 = np.isin(event.pid, decay_handler.all_stable_pids) & (
            event.status == 1
        )
        # Assert that number of stable particles of given type hasn't decreased
        assert np.sum(stable0) <= np.sum(stable1)
        # Assert that stable particles haven't decayed
        assert np.all(stable0 == stable1[0 : len(event0)])

        # Assert that initial stack is preserved in decayed stack
        # Fix after decay
        event_dec.status = event0.status
        event_dec.daughters = event0.daughters
        # event_dec is np.float32: 1/3, 2/3 are represented differently
        # from event0.charge (np.float64), therefore:
        event0.charge = np.array(event0.charge, dtype=np.float32)

        # Fix after selection
        event_dec.mothers = event.mothers[0 : len(event0)]
        assert event_dec == event0


@pytest.mark.parametrize("Model", get_all_models())
def test_decay_handler(Model):
    stable_particles = [111]
    decaying_particles = [-211, 211, 13, -13, 2112]

    if Model.name in ["UrQMD", "QGSJet"]:
        # UrQMD produce neutrons = 2112 with wrong masses
        # QGSJet produce neutrons with proton mass
        # Pythia8 ignores them for decay
        decaying_particles = [-211, 211, 13, -13]

    if Model.name in ["PhoJet", "Pythia"]:
        evt_kin = chromo.kinematics.FixedTarget(100, "p", "p")

    elif Model.name == "Sophia":
        evt_kin = chromo.kinematics.FixedTarget(100, "gamma", "p")
    else:
        evt_kin = chromo.kinematics.FixedTarget(100, "p", "O")

    run_in_separate_process(
        run_decay_handler, Model, evt_kin, stable_particles, decaying_particles
    )
