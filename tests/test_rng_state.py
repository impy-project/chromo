import impy
from impy.kinematics import FixedTarget, CenterOfMass
from impy.constants import TeV, GeV
import impy.models as im
import pickle
from pathlib import Path

import abc
import pytest
from .util import run_in_separate_process, xfail_on_ci_if_model_is_incompatible

# generate list of all models in impy.models
models = list(obj for obj in im.__dict__.values() if type(obj) is abc.ABCMeta)


def rng_state_test(model):

    if model is im.Sophia20:
        ekin = FixedTarget(13 * TeV, "photon", "proton")
    elif model is im.UrQMD34:
        ekin = CenterOfMass(50 * GeV, "proton", "proton")
    else:
        ekin = CenterOfMass(13 * TeV, "proton", "proton")

    generator = model(ekin, seed=3163325)
    nevents = 10
    rng_state_file = str(model.name) + "rng_state.dat"

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
        assert (
            counters[i] == counter
        ), 'Counters for seed {0} after "{1}" event are different:\n expected(previous) = {2}, received(current) = {3}'.format(
            generator.random_state.seed, i, counters[i], counter
        )
        i = i + 1

    # Test for restoring state from file:
    state_after_now = generator.random_state.copy()
    # Restore from file
    with open(rng_state_file, "rb") as pfile:
        generator.random_state = pickle.load(pfile)

    if Path(rng_state_file).exists():
        Path(rng_state_file).unlink()

    # And check for equality
    state_equal = None
    if state_after_now == generator.random_state:
        state_equal = True
    else:
        state_equal = False
    assert state_equal, "Restored state from file is different from obtained"


@pytest.mark.parametrize("Model", models)
def test_generator(Model):
    xfail_on_ci_if_model_is_incompatible(Model)
    run_in_separate_process(rng_state_test, Model)
