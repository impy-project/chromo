from pathlib import Path
from impy.constants import GeV
from impy.kinematics import CenterOfMass, FixedTarget
import impy.models as im
import pytest
import pyhepmc
import abc
from .util import run_in_separate_process, xfail_on_ci_if_model_is_incompatible

# generate list of all models in impy.models
models = set(obj for obj in im.__dict__.values() if type(obj) is abc.ABCMeta)


def run(Model):
    ekin = CenterOfMass(10 * GeV, 2212, 2212)
    if Model == im.Sophia20:
        ekin = FixedTarget(10 * GeV, "photon", "proton")
    gen = Model(ekin, seed=715434)
    return list(ev.copy() for ev in gen(3))


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
    print(test_file)

    xfail_on_ci_if_model_is_incompatible(Model)
    if Model in (
        im.DpmjetIII306,
        im.DpmjetIII191,
        im.DpmjetIII193,
        im.Phojet112,
        im.Phojet191,
        im.UrQMD34,
    ):
        pytest.xfail("The model fails for some parameters")

    events = run_in_separate_process(run, Model, timeout=60)

    expected = []
    with pyhepmc.open(test_file, "w") as f:
        for event in events:
            f.write(event)
            expected.append(event.to_hepmc3())

    restored = []
    with pyhepmc.open(test_file) as events:
        for event in events:
            assert event is not None
            # Some models (e.g. Sibyll) do not set event.event_number properly
            # assert event.event_number == ievent
            restored.append(event)

    assert len(restored) == len(expected)

    for ev1, ev2 in zip(expected, restored):
        assert ev1 == ev2

    # delete test_file if test is successful
    test_file.unlink()
