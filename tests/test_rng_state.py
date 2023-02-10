from chromo.kinematics import FixedTarget, CenterOfMass
from chromo.constants import TeV, GeV, PeV
import chromo.models as im
import pickle
import pytest
from chromo.util import get_all_models, name2pdg, CompositeTarget
from copy import deepcopy
import numpy as np


proton = name2pdg("proton")
photon = name2pdg("photon")
composite = CompositeTarget([("p", 0.5), ("N", 0.5)])


@pytest.mark.parametrize("Model", get_all_models())
@pytest.mark.parametrize("target", ("proton", "composite"))
def test_rng_state_1(Model, target):
    if Model is im.Pythia8:
        pytest.skip("Pythia8 currently does not support rng_state serialization")

    if Model is im.UrQMD34:
        pytest.xfail("UrQMD fails this test for most seeds, needs investigation")

    p = proton if proton in Model.projectiles else photon
    t = globals()[target]

    if t not in Model.targets:
        pytest.skip(f"{Model.pyname} does not support target {target}")

    kin = FixedTarget(1 * PeV, p, t)

    generator = Model(seed=1)
    nevents = 10

    # Save a initial state to a variable:
    state_0 = deepcopy(generator.random_state)

    # Generate nevents events
    states = []
    for _ in generator(kin, nevents):
        states.append(deepcopy(generator.random_state))

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
@pytest.mark.parametrize("new_model", (False, True))
def test_rng_state_2(Model, new_model):
    if Model is im.Pythia8:
        pytest.skip("Pythia8 currently does not support rng_state serialization")

    if Model in (im.UrQMD34, im.EposLHC):
        pytest.xfail(f"{Model.pyname} fails this test, needs investigation")

    if Model is im.Sophia20:
        kin = FixedTarget(13 * TeV, "photon", "proton")
    elif Model is im.UrQMD34:
        kin = CenterOfMass(50 * GeV, "proton", "proton")
    else:
        kin = CenterOfMass(13 * TeV, "proton", "proton")

    model = Model(seed=1)

    # Save a initial state to a variable:
    state_0 = model.random_state.copy()

    # Generate some events
    events1 = []
    nevents = 10
    for event in model(kin, nevents):
        events1.append(event.copy())

    # Pickle generator state after nevents
    pickled_state_1 = pickle.dumps(model.random_state)

    if new_model:
        del model
        # seed should be irrelevant when we restore state
        model = Model(seed=12345)

    # Restore initial state from variable
    model.random_state = state_0

    # Generate some more events
    events2 = []
    for event in model(kin, nevents):
        events2.append(event.copy())

    state_2 = model.random_state
    state_1 = pickle.loads(pickled_state_1)

    # check that pickled state_1 and reproduced state_2 are equal
    assert state_2 == state_1

    # check that sequences of generated events are equal
    # comparing pid arrays is enough
    events1_pid = np.concatenate([x.pid for x in events1])
    events2_pid = np.concatenate([x.pid for x in events2])

    np.testing.assert_equal(events1_pid, events2_pid)
