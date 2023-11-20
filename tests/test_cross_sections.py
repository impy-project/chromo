from chromo.constants import GeV
from chromo.kinematics import (
    CenterOfMass,
    CompositeTarget,
)
import chromo.models as im
import pytest
from .util import run_in_separate_process
from chromo.util import get_all_models
import dataclasses
import numpy as np


def run_model(Model, kin, number=1):
    try:
        gen = Model(kin, seed=1)
    except ValueError:
        return None

    return gen.cross_section()


@pytest.mark.parametrize("target", ("p", "air"))
@pytest.mark.parametrize("projectile", ("gamma", "pi-", "K+", "p", "He"))
@pytest.mark.parametrize("Model", get_all_models())
def test_generator(projectile, target, Model):
    p1 = projectile
    p2 = target
    if Model is im.UrQMD34:
        pytest.skip("UrQMD34 cross sections broken.")
    if p2 == "air":
        if Model is im.Pythia8:
            pytest.skip("Simulating nuclei in Pythia8 is very time-consuming")

        # cannot use Argon in SIBYLL, so make air from N, O only
        p2 = CompositeTarget((("N", 0.78), ("O", 0.22)))

    kin = CenterOfMass(100 * GeV, p1, p2)

    cs = run_in_separate_process(run_model, Model, kin)
    if cs is None:
        assert abs(kin.p1) not in Model.projectiles or abs(kin.p2) not in Model.targets
        return
    cstup = dataclasses.astuple(cs)
    assert np.any(cstup), "No valid cross sections returned"
    if p1 == "gamma":
        assert np.nansum(cstup) > 0.1, "Cross section too small" + str(cs)
    else:
        assert np.nansum(cstup) > 20, "Cross section too small" + str(cs)
