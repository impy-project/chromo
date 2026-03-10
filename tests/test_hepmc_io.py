from contextlib import nullcontext
from pathlib import Path

import pyhepmc
import pytest
from pytest import approx

import chromo.models as im
from chromo.constants import GeV
from chromo.kinematics import CenterOfMass, FixedTarget
from chromo.util import get_all_models

from .util import run_in_separate_process

# generate list of all models in chromo.models
models = get_all_models()


def run(Model):
    evt_kin = CenterOfMass(100 * GeV, "proton", "proton")

    if Model == im.Sophia20:
        evt_kin = FixedTarget(100 * GeV, "photon", "proton")
    gen = Model(evt_kin, seed=1)
    inel_cs = gen.cross_section().inelastic
    return list(gen(3)), inel_cs


@pytest.mark.parametrize(
    "Model",
    models,
)
def test_hepmc_io(Model):
    # To run this test do `pytest tests/test_hepmc_writer.py`
    # This test fails because the event record written by HepMC3 C++ is bad,
    # a lot of particles are missing. Either a bug in the original chromo record or a
    # bug in the HepMC3 C++ code (not the pyhepmc code).

    test_file = Path(f"{Path(__file__).stem}_{Model.pyname}.dat")

    events, inel_cs = run_in_separate_process(run, Model)
    expected = []
    genevent = None

    with (
        pytest.warns(RuntimeWarning)
        if (Model.name in ["DPMJET-III", "PhoJet"] or Model.pyname == "Pythia8")
        else nullcontext()
    ):
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
            restored.append(event)  # noqa: PERF402

    assert len(restored) == len(expected)

    for ev1, ev2 in zip(expected, restored):
        assert ev1 == ev2
        if Model.name.startswith("UrQMD"):
            continue
        assert ev1.cross_section.xsec() / 1e9 == approx(inel_cs)
        assert ev2.cross_section.xsec() / 1e9 == approx(inel_cs)

    # delete test_file if test is successful
    test_file.unlink()
