from pathlib import Path
from impy.constants import GeV
from impy.kinematics import CenterOfMass, FixedTarget
import impy.models as im
import pytest
import pyhepmc
from impy.util import get_all_models

# generate list of all models in impy.models
models = get_all_models()


@pytest.mark.parametrize(
    "Model",
    models,
)
def test_hepmc_io(Model):
    # To run this test do `pytest tests/test_hepmc_writer.py`
    # This test fails because the event record written by HepMC3 C++ is bad,
    # a lot of particles are missing. Either a bug in the original impy record or a
    # bug in the HepMC3 C++ code (not the pyhepmc code).

    test_file = Path(f"{Path(__file__).stem}_{Model.pyname}.dat")

    kin = CenterOfMass(100 * GeV, "proton", "proton")

    if Model == im.Sophia20:
        kin = FixedTarget(100 * GeV, "photon", "proton")
    gen = Model(seed=1)
    events = [x.copy() for x in gen(kin, 3)]

    expected = []
    genevent = None
    for ev in events:
        genevent = ev.to_hepmc3(genevent)
        expected.append(genevent)

    # Uncomment this to get debugging output. Higher number shows more.
    # This only works if you compile pyhepmc in debug mode.
    # pyhepmc.Setup.debug_level = 100

    with pyhepmc.open(test_file, "w") as f:
        for event in expected:
            f.write(event)

    restored = []
    with pyhepmc.open(test_file) as f:
        for event in f:
            restored.append(event)

    assert len(restored) == len(expected)

    for ev1, ev2 in zip(expected, restored):
        assert ev1 == ev2

    # delete test_file if test is successful
    test_file.unlink()
