from impy.constants import TeV
from impy.kinematics import EventKinematics
import impy.models as im
import abc
from collections import Counter
import pytest
from multiprocessing import Pool
from multiprocessing.context import TimeoutError
import os
from .util import run_in_separate_process

# generate list of all models in impy.models
models = set(obj for obj in im.__dict__.values() if type(obj) is abc.ABCMeta)


def run_model(model, ekin):
    gen = model(ekin)

    c = Counter()
    for event in gen(10):
        ev = event.final_state()
        assert len(ev.pid) > 0
        c.update(ev.pid)

    return c


@pytest.mark.parametrize("Model", models)
def test_generator(Model):
    # remove this when git lfs issue is fixed
    if os.environ.get("CI", False) and Model in (
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
    if Model is im.Sophia20:
        # Sophia can only do γp, γn
        p1pdg = 22  # gamma
    elif Model in [im.Phojet112, im.UrQMD34]:
        # The old phojet needs more tweaking for pion-proton (is not related to test)
        p1pdg = 2212  # proton

    ekin = EventKinematics(
        ecm=7 * TeV,
        p1pdg=p1pdg,
        p2pdg=p2pdg,
    )

    c = run_in_separate_process(run_model, Model, ekin)

    assert c[211] > 0, "pi+"
    assert c[-211] > 0, "pi-"
    assert c[2212] > 0, "p"
