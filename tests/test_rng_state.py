from impy.kinematics import FixedTarget, CenterOfMass
from impy.constants import TeV, GeV
import impy.models as im
import pickle
from pathlib import Path
import pytest
from .util import (
    run_in_separate_process,
    get_all_models,
)

# generate list of all models in impy.models
Models = get_all_models(im)


def run_rng_state(Model):
    if Model is im.Sophia20:
        evt_kin = FixedTarget(13 * TeV, "photon", "proton")
    elif Model is im.UrQMD34:
        evt_kin = CenterOfMass(50 * GeV, "proton", "proton")
    else:
        evt_kin = CenterOfMass(13 * TeV, "proton", "proton")

    generator = Model(evt_kin, seed=1)
    nevents = 10
    rng_state_file = Path(f"{Model.pyname}_rng_state.dat")

    # Save a initial state to a variable:
    state0 = generator.random_state.copy()

    # Generate nevents events
    counters = []
    for event in generator(nevents):
        counters.append(generator.random_state.counter)

    # Save generator state after nevents to a file
    with open(rng_state_file, "wb") as pfile:
        pickle.dump(generator.random_state, pfile, protocol=pickle.HIGHEST_PROTOCOL)

    # Restore initial state from variable
    generator.random_state = state0

    # And compare counters after each generated event
    i = 0
    for event in generator(nevents):
        counter = generator.random_state.counter
        assert counters[i] == counter, (
            f"Counters for seed {generator.random_state.seed} "
            f"after {i} event are different\n"
            f"  expected(previous) = {counters[i]}, received(current) = {counter}"
        )
        i = i + 1

    # Test for restoring state from file:
    state_after_now = generator.random_state.copy()
    # Restore from file
    with open(rng_state_file, "rb") as pfile:
        generator.random_state = pickle.load(pfile)

    rng_state_file.unlink()

    # And check for equality
    state_equal = None
    if state_after_now == generator.random_state:
        state_equal = True
    else:
        state_equal = False
    assert state_equal, "Restored state from file is different from obtained"


@pytest.mark.parametrize("Model", Models)
def test_rng_state(Model):
    if Model is im.Pythia8:
        pytest.skip("Pythia8 currently does not support rng_state serialization")
    run_in_separate_process(run_rng_state, Model)
