import chromo
import numpy as np
import pytest
import sys
import os
import platform

from chromo.constants import long_lived
from chromo.decay_handler import Pythia8DecayHandler
from chromo.util import get_all_models
from .util import run_in_separate_process

if sys.platform == "win32":
    pytest.skip("Pythia8 does not run on windows", allow_module_level=True)


def run_decay_handler(model, evt_kin, stable_particles):
    seed = 1
    decay_handler = Pythia8DecayHandler(stable_particles, seed)
    gen = model(evt_kin, seed=seed)

    for event in gen(10):
        event0 = event.copy()
        decay_handler(event)

        # Assert that decay has happend
        event_dec = event[0 : len(event0)]
        assert np.sum(event_dec.status != 1) >= np.sum(event0.status != 1)

        # Assert that all particles that should decay have been decayed
        not_decayed0 = np.isin(event0.pid, decay_handler.all_unstable_pids) & (
            event0.status == 1
        )
        not_decayed1 = np.isin(event.pid, decay_handler.all_unstable_pids) & (
            event.status == 1
        )

        if np.sum(not_decayed0) > 0:
            failed_to_decay = event.pid[not_decayed1]
            assert (
                np.sum(not_decayed1) == 0
            ), f"{failed_to_decay} with status {event.status[not_decayed1]} do not decay"

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

        # Assert that initial stack is preserved in decayed stack, i.e.
        # original particles haven't been changed, except for
        # status and daughters

        # We want to use comparison operation which compares
        # all fields. Because some fields have changed, assign them
        # to original values to pass assertion
        event_dec.status = event0.status
        event_dec.daughters = event0.daughters
        # event_dec is np.float32: 1/3, 2/3 are represented differently
        # from event0.charge (np.float64), therefore:
        event0.charge = np.array(event0.charge, dtype=np.float32)

        # `event_dec = event[0 : len(event0)]` operation changes
        # mothers by using `_select_mothers`
        # revert this so that comparison
        # `event_dec == event0` will pass assertion
        event_dec.mothers = event.mothers[0 : len(event0)]
        assert event_dec == event0


@pytest.mark.parametrize("Model", get_all_models())
def test_decay_handler(Model):
    stable_particles = long_lived

    if (
        (os.environ.get("CI", False))
        and (platform.system() == "Darwin")
        and (Model.name == "UrQMD")
    ):
        pytest.xfail(
            f"For {Model.pyname} DecayHandler fails to decay all unstable particles on MacOS CI"
        )

    if Model.name in ["PhoJet", "Pythia"]:
        evt_kin = chromo.kinematics.FixedTarget(100, "p", "p")

    elif Model.name == "Sophia":
        evt_kin = chromo.kinematics.FixedTarget(100, "gamma", "p")
    else:
        evt_kin = chromo.kinematics.FixedTarget(100, "p", "O")

    run_in_separate_process(run_decay_handler, Model, evt_kin, stable_particles)
