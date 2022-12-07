from impy.kinematics import FixedTarget, CenterOfMass
from impy.constants import TeV, GeV
import impy.models as im
import pickle
import pytest
from .util import run_in_separate_process
from impy.util import get_all_models

# generate list of all models in impy.models
Models = get_all_models()


def run_rng_state(Model):
    if Model is im.Sophia20:
        evt_kin = FixedTarget(13 * TeV, "photon", "proton")
    elif Model is im.UrQMD34:
        evt_kin = CenterOfMass(50 * GeV, "proton", "proton")
    else:
        evt_kin = CenterOfMass(13 * TeV, "proton", "proton")

    generator = Model(evt_kin, seed=1)
    nevents = 10

    # Save a initial state to a variable:
    state_0 = generator.random_state.copy()

    # Generate nevents events
    counters = []
    for _ in generator(nevents):
        counters.append(generator.random_state.counter)

    # Save generator state after nevents
    pickled_state_1 = pickle.dumps(generator.random_state)

    # Restore initial state from variable
    generator.random_state = state_0

    # And compare counters after each generated event
    i = 0
    for _ in generator(nevents):
        counter = generator.random_state.counter
        assert counters[i] == counter, (
            f"Counters for seed {generator.random_state.seed} "
            f"after {i} event are different\n"
            f"  expected(previous) = {counters[i]}, received(current) = {counter}"
        )
        i = i + 1

    state_2 = generator.random_state.copy()
    state_1 = pickle.loads(pickled_state_1)

    # check pickled state_1 and reproduced state_2 for equality
    assert state_2 == state_1, "Restored state from file is different from obtained"


@pytest.mark.parametrize("Model", Models)
def test_rng_state(Model):
    if Model is im.Pythia8:
        pytest.skip("Pythia8 currently does not support rng_state serialization")
    run_in_separate_process(run_rng_state, Model)
