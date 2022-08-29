from impy.constants import TeV
from impy.kinematics import EventKinematics
import impy.models as im
import abc
from collections import Counter
import pytest
from multiprocessing import Pool
from multiprocessing.context import TimeoutError
import os

# generate list of all models in impy.models
models = set(obj for obj in im.__dict__.values() if type(obj) is abc.ABCMeta)

# models = [im.DpmjetIII306, im.Phojet112]

def run_model(model, ekin):
    gen = model(ekin)

    c = Counter()
    for event in gen(10):
        ev = event.final_state()
        assert len(ev.pid) > 0
        c.update(ev.pid)

    return c


# DpmjetIII191 sometimes randomly fails this
@pytest.mark.parametrize("model", models)
def test_generators(model):
    # remove this when git lfs issue is fixed
    if os.environ.get("CI", False) and model in (
        im.QGSJet01c,
        im.QGSJetII03,
        im.QGSJetII04,
        im.Phojet191,
        im.EposLHC,
        im.DpmjetIII306,
        im.DpmjetIII191,
        im.DpmjetIII193,
    ):
        pytest.xfail("model cannot succeed on CI, because git lfs does not work")

    p1pdg = -211  # pi-
    p2pdg = 2212  # proton
    if model is im.Sophia20:
        # Sophia can only do γp, γn
        p1pdg = 22  # gamma
    elif model in [im.Phojet112, im.UrQMD34]:
        # The old phojet needs more tweaking for pion-proton (is not related to test)
        p1pdg = 2212  # proton

    ekin = EventKinematics(
        ecm=7 * TeV,
        p1pdg=p1pdg,
        p2pdg=p2pdg,
    )

    # Some models need to initialize same fortran code, which can only be
    # initialized once. As a workaround, we run each model in a separate
    # thread. When running several jobs, maxtasksperchild=1 is needed to
    # use a fresh interpreter for each task (not needed here, but still).
    with Pool(1, maxtasksperchild=1) as p:
        r = p.apply_async(run_model, (model, ekin))
        try:
            c = r.get(timeout=30)
        except TimeoutError:
            # usually happens when model aborts and kills child process
            raise TimeoutError("check stdout for errors")

    assert c[211] > 0, "pi+"
    assert c[-211] > 0, "pi-"
    assert c[2212] > 0, "p"
