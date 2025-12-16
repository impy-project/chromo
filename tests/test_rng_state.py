import pickle
from copy import deepcopy

import numpy as np
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


def run_rng_state_with_bitgen(Model, bitgen_class, seed):
    """Test RNG state save/restore with specific bit generator."""
    if Model is im.Sophia20:
        evt_kin = FixedTarget(13 * TeV, "photon", "proton")
    elif Model is im.UrQMD34:
        evt_kin = CenterOfMass(50 * GeV, "proton", "proton")
    else:
        evt_kin = CenterOfMass(13 * TeV, "proton", "proton")

    rng = np.random.Generator(bitgen_class(seed=seed))
    generator = Model(evt_kin, seed=rng)

    nevents = 10 if Model is not im.EposLHC else 50

    state_0 = deepcopy(generator.random_state)
    states = [deepcopy(generator.random_state) for _ in generator(nevents)]

    generator.random_state = state_0

    for i, _ in enumerate(generator(nevents)):
        current_state = generator.random_state
        # Use str() for comparison to handle numpy arrays in nested dicts
        assert str(states[i]) == str(
            current_state
        ), f"states differ after {i+1} generation with {bitgen_class.__name__}"


@pytest.mark.parametrize("Model", get_all_models())
def test_rng_state(Model):
    if Model in (im.UrQMD34,):
        #       UrQMD has internal state that affects event
        #       generation but is not captured by RNG state
        #       There are several places where it can occur:
        #        - Commong blocks: I tried to save and restore all commong blocks
        #        with help of Copilot. It doesn't work, meaning that answer buried deep
        #       - Save statements: they are not reset and might influence rng
        #       - Pythia internal rng: it might be but, it seems pythia uses correct rng
        #       Noticed: there are kind of "warm-up" that calculates some variables
        #       that at later stages is not called. There are many different rng in the code
        #       it might be that connecting all of them to one rng merge different rng branches
        #       that must be independent: event and variable sequencies
        pytest.xfail(f"{Model.pyname} fails this test, needs investigation")
    if Model in (im.EposLHCR, im.EposLHCRHadrRescattering):
        pytest.skip(
            f"{Model.pyname} maintains UrQMD internal state that affects event "
            "generation but is not captured by RNG state."
        )

    run_in_separate_process(run_rng_state, Model)


@pytest.mark.parametrize("Model", get_all_models())
@pytest.mark.parametrize(
    "bitgen_class,seed",
    [
        (np.random.PCG64, 11111),
        (np.random.MT19937, 22222),
        (np.random.Philox, 33333),
        (np.random.SFC64, 44444),
        (np.random.PCG64DXSM, 55555),
    ],
)
def test_rng_state_bitgens(Model, bitgen_class, seed):
    """Test different NumPy bit generators with all models."""
    if Model in (im.UrQMD34,):
        pytest.xfail(f"{Model.pyname} fails this test, needs investigation")
    if Model in (im.EposLHCR, im.EposLHCRHadrRescattering):
        pytest.skip(
            f"{Model.pyname} maintains UrQMD internal state that affects event "
            "generation but is not captured by RNG state."
        )

    # Skip Pythia6 + MT19937 combination - model-specific interaction issue
    if Model is im.Pythia6 and bitgen_class == np.random.MT19937:
        pytest.skip(
            "Pythia6 has unknown issues with MT19937 bit generator, needs investigation"
        )

    run_in_separate_process(run_rng_state_with_bitgen, Model, bitgen_class, seed)
