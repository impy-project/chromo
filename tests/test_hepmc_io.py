from pathlib import Path
from impy.constants import GeV
from impy.kinematics import CenterOfMass, FixedTarget
import impy.models as im
import pytest
import pyhepmc
from .util import (
    run_in_separate_process,
    xfail_on_ci_if_model_is_incompatible,
    get_all_models,
)

# generate list of all models in impy.models
models = get_all_models(im)


def run(Model):
    ekin = CenterOfMass(10 * GeV, "proton", "proton")
    if Model == im.Sophia20:
        ekin = FixedTarget(10 * GeV, "photon", "proton")
    gen = Model(ekin, seed=1)
    return list(gen(3))


@pytest.mark.parametrize(
    "Model",
    models,
)
def test_hepmc_io(Model):
    # To run this test do `pytest tests/test_hepmc_writer.py`
    # This test fails because the event record written by HepMC3 C++ is bad,
    # a lot of particles are missing. Either a bug in the original impy record or a
    # bug in the HepMC3 C++ code (not the pyhepmc code).

    test_file = Path(f"{Path(__file__).with_suffix('')}_{Model.__name__}.dat")

    xfail_on_ci_if_model_is_incompatible(Model)

    events = run_in_separate_process(run, Model, timeout=30)
    expected = [ev.to_hepmc3() for ev in events]

    with pyhepmc.open(test_file, "w") as f:
        for event in expected:
            f.write(event)  # this calls to_hepmc3() implicitly

    restored = []
    with pyhepmc.open(test_file) as f:
        for event in f:
            restored.append(event)

    assert len(restored) == len(expected)

    for ev1, ev2 in zip(expected, restored):
        assert ev1 == ev2

    # delete test_file if test is successful
    test_file.unlink()
