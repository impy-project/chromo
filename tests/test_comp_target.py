from impy.constants import TeV, GeV
from impy.kinematics import CenterOfMass, CompositeTarget
import impy.models as im
from collections import Counter
import pytest
from multiprocessing import Pool
from multiprocessing.context import TimeoutError
from .util import xfail_on_ci_if_model_is_incompatible, get_all_models

# generate list of all models in impy.models
models = get_all_models(im)


def run_model(model, ekin):
    gen = model(ekin, seed=1)

    c = Counter()
    for event in gen(10):
        ev = event.final_state()
        assert len(ev.pid) > 0
        c.update(ev.pid)

    return c


@pytest.mark.parametrize("model", models)
def test_generators(model):
    xfail_on_ci_if_model_is_incompatible(model)
    if model in (
        im.Sophia20,
        im.Phojet112,
        im.Phojet191,
        im.Phojet193,
    ):
        pytest.xfail("Model doesn't support nuclei")

    if model in (im.Sibyll23d,):
        pytest.xfail("Sibyll23d is brocken. Something with rng.")

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
    ekin = CenterOfMass(7 * TeV, projectile, target)

    if model in (im.UrQMD34,):
        ekin = CenterOfMass(50 * GeV, projectile, target)

    # Some models need to initialize same fortran code, which can only be
    # initialized once. As a workaround, we run each model in a separate
    # thread. When running several jobs, maxtasksperchild=1 is needed to
    # use a fresh interpreter for each task (not needed here, but still).
    with Pool(1, maxtasksperchild=1) as p:
        r = p.apply_async(run_model, (model, ekin))
        try:
            c = r.get(timeout=60)
        except TimeoutError:
            # usually happens when model aborts and kills child process
            raise TimeoutError("check stdout for errors")

    assert c[211] > 0, "pi+"
    assert c[-211] > 0, "pi-"
    assert c[2212] > 0, "p"
