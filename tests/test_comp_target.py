from impy.constants import TeV, GeV
from impy.kinematics import CenterOfMass, CompositeTarget
import impy.models as im
from collections import Counter
import pytest
from multiprocessing import Pool
from multiprocessing.context import TimeoutError
from .util import (
    xfail_on_ci_if_model_is_incompatible,
    get_all_models,
    run_in_separate_process,
)

# generate list of models to test,
# skip models which do not support nuclei
Models = get_all_models(im)
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
    xfail_on_ci_if_model_is_incompatible(Model)

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
    evt_kin = CenterOfMass(7 * TeV, projectile, target)

    if model in (im.UrQMD34,):
        evt_kin = CenterOfMass(50 * GeV, projectile, target)

    # Some models need to initialize same fortran code, which can only be
    # initialized once. As a workaround, we run each model in a separate
    # thread. When running several jobs, maxtasksperchild=1 is needed to
    # use a fresh interpreter for each task (not needed here, but still).
    with Pool(1, maxtasksperchild=1) as p:
        r = p.apply_async(run_model, (model, evt_kin))
        try:
            c = r.get(timeout=60)
        except TimeoutError:
            # usually happens when model aborts and kills child process
            raise TimeoutError("check stdout for errors")

    assert c[211] > 0, "pi+"
    assert c[-211] > 0, "pi-"
    assert c[2212] > 0, "p"
