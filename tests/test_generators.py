from impy.constants import GeV
from impy.kinematics import CenterOfMass
import impy.models as im
from collections import Counter
import pytest
from .util import run_in_separate_process
from impy.util import get_all_models


# generate list of all models in impy.models
models = get_all_models()


def run_model(Model, evt_kin):
    gen = Model(evt_kin, seed=1)

    c = Counter()
    for event in gen(10):
        ev = event.final_state()
        assert len(ev.pid) > 0
        c.update(ev.pid)

    return c


@pytest.mark.parametrize("Model", models)
def test_generator(Model):
    p1pdg = -211  # pi-
    p2pdg = 2212  # proton
    if Model is im.Sophia20:
        # Sophia can only do γp, γn
        p1pdg = 22  # gamma
    elif Model in [im.Phojet112, im.UrQMD34]:
        # The old phojet needs more tweaking for pion-proton (is not related to test)
        p1pdg = 2212  # proton

    evt_kin = CenterOfMass(
        100 * GeV,
        p1pdg,
        p2pdg,
    )

    c = run_in_separate_process(run_model, Model, evt_kin)

    assert c[211] > 0, "pi+"
    assert c[-211] > 0, "pi-"
    assert c[2212] > 0, "p"
