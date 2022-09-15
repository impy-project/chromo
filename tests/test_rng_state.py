import impy
from impy.kinematics import FixedTarget, CenterOfMass
from impy.constants import TeV, GeV
import impy.models as im

import abc
import pytest
from .util import run_in_separate_process, xfail_on_ci_if_model_is_incompatible
import os

# generate list of all models in impy.models
models = set(obj for obj in im.__dict__.values() if type(obj) is abc.ABCMeta)


def rng_state_test(model):

    if model is im.Sophia20:
        ekin = FixedTarget(13 * TeV, "photon", "proton")
    elif model is im.UrQMD34:
        ekin = CenterOfMass(50 * GeV, "proton", "proton")
    else:
        ekin = CenterOfMass(13 * TeV, "proton", "proton")

    generator = model(ekin, seed=3163325)
    nevents = 10

    # Save a initial state to a variable:
    state0 = generator.rng_state.copy()

    # Generate nevents events
    counters = []
    for event in generator(nevents):
        counters.append(generator.rng_state.counter)

    # Save generator state after nevents to a file
    generator.dump_rng_state_to("rng_state.dat")

    # Restore initial state from variable
    generator.rng_state = state0

    # And compare counters after each generated event
    i = 0
    for event in generator(nevents):
        counter = generator.rng_state.counter
        assert (
            counters[i] == counter
        ), 'Counters after "{0}" event are different:\n expected(previous) = {1}, received(current) = {2}'.format(
            i, counters[i], counter
        )
        i = i + 1

    # Test for restoring state from file:
    state_after_now = generator.rng_state.copy()
    # Restore from file
    generator.restore_rng_state_from("rng_state.dat")

    if os.path.exists("rng_state.dat"):
        os.remove("rng_state.dat")

    # And check for equality
    state_equal = None
    if state_after_now == generator.rng_state:
        state_equal = True
    else:
        state_equal = False
    assert state_equal, "Restored state from file is different from obtained"


@pytest.mark.parametrize("Model", models)
def test_generator(Model):
    xfail_on_ci_if_model_is_incompatible(Model)
    run_in_separate_process(rng_state_test, Model)


# rng_state_test(im.UrQMD34)
