import pickle
from copy import deepcopy

import pytest

import chromo.models as im
from chromo.constants import GeV, TeV
from chromo.kinematics import CenterOfMass, FixedTarget
from chromo.util import get_all_models

from .util import run_in_separate_process


def run_rng_state(Model):
    if Model is im.Sophia20:
        evt_kin = FixedTarget(13 * TeV, "photon", "proton")
    elif Model is im.UrQMD34:
        evt_kin = CenterOfMass(50 * GeV, "proton", "proton")
    else:
        evt_kin = CenterOfMass(13 * TeV, "proton", "proton")

    generator = Model(evt_kin, seed=1)

    nevents = 10
    # EposLHC needs more events to fail
    if Model is im.EposLHC:
        nevents = 50

    # Save a initial state to a variable:
    state_0 = deepcopy(generator.random_state)
    # Generate nevents events
    states = [deepcopy(generator.random_state) for _ in generator(nevents)]

    # Pickle generator state after nevents
    pickled_state_1 = pickle.dumps(generator.random_state)

    # Restore initial state from variable
    generator.random_state = state_0

    # Compare states after each generated event
    for i, _ in enumerate(generator(nevents)):
        assert (
            states[i] == generator.random_state
        ), f"states differ after {i+1} generation"

    state_1a = generator.random_state
    state_1 = pickle.loads(pickled_state_1)
    # check that pickled state_1 and reproduced state_1a are equal
    assert state_1a == state_1


@pytest.mark.parametrize("Model", get_all_models())
def test_rng_state(Model):
    if Model is im.Pythia8:
        pytest.skip("Pythia8 currently does not support rng_state serialization")

    if Model in (im.UrQMD34,):
        pytest.xfail(f"{Model.pyname} fails this test, needs investigation")

    run_in_separate_process(run_rng_state, Model)
