import dataclasses

import numpy as np
import pytest

import chromo.models as im
from chromo.constants import GeV
from chromo.kinematics import (
    CenterOfMass,
    CompositeTarget,
)
from chromo.util import get_all_models

from .util import run_in_separate_process


def run_model(Model, kin):
    try:
        gen = Model(kin, seed=1)
    except ValueError:
        return None

    return (gen.cross_section(), gen._inel_or_prod_cross_section)


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

    # increase energy for K+ and SIBYLL23c to avoid prod cs = nan region
    # this is specifig to SIBYLL23C
    if Model == im.Sibyll23c and projectile == "K+" and target == "air":
        kin = CenterOfMass(1e3 * GeV, p1, p2)
    else:
        kin = CenterOfMass(100 * GeV, p1, p2)

    ret = run_in_separate_process(run_model, Model, kin)
    if ret is None:
        assert abs(kin.p1) not in Model.projectiles or abs(kin.p2) not in Model.targets
        return
    cs, private_prod = ret
    cstup = dataclasses.astuple(cs)
    assert np.any(cstup), "No valid cross sections returned"
    if p1 == "gamma":
        assert np.nansum(cstup) > 0.1, "Cross section too small" + str(cs)
    else:
        assert np.nansum(cstup) > 20, "Cross section too small" + str(cs)

    assert private_prod is not None and private_prod > 0, "No private cross section set"
