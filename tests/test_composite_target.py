from impy.constants import TeV, GeV
from impy.kinematics import CenterOfMass, CompositeTarget
import impy.models as im
from collections import Counter
import pytest
from multiprocessing import Pool
from multiprocessing.context import TimeoutError
from .util import (
    skip_on_ci_if_model_is_incompatible,
    get_all_models,
    run_in_separate_process,
)

# generate list of models to test,
Models = get_all_models(im)
# skip models which do not support nuclei
Models = [M for M in Models if M.name not in ("Sophia", "PhoJet")]


def run_model(Model, evt_kin):
    gen = Model(evt_kin, seed=1)

    c = Counter()
    for event in gen(10):
        ev = event.final_state()
        assert len(ev.pid) > 0
        c.update(ev.pid)

    return c


@pytest.mark.parametrize("Model", Models)
def test_composite_target(Model):
    skip_on_ci_if_model_is_incompatible(Model)

    if Model is im.Pythia8:
        pytest.skip("Switching beams in Pythia8 is very time-consuming")

    projectile = "pi-"
    seed_for_test = 321
    target = CompositeTarget(
        [
            ("N14", 2 * 0.78084, "Nitrogen"),
            ("O16", 2 * 0.20946, "Oxygen"),
            ("O16", 0.004, "Oxygen(Vapor)"),
            ("proton", 2 * 0.004, "Hydrogen(Vapor)"),
        ],
        "Air without argon",
        seed_for_test,
    )

    evt_kin = CenterOfMass(10 * GeV, projectile, target)

    c = run_in_separate_process(run_model, Model, evt_kin)

    assert c[211] > 0, "pi+"
    assert c[-211] > 0, "pi-"
    assert c[2212] > 0, "p"
