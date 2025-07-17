import numpy as np
import pytest
from numpy.testing import assert_allclose

import chromo.models as im
from chromo.constants import GeV
from chromo.kinematics import (
    CenterOfMass,
    CompositeTarget,
    EventFrame,
    EventKinematicsMassless,
    EventKinematicsWithRestframe,
    FixedTarget,
    TotalEnergy,
)
from chromo.util import get_all_models

from .util import run_in_separate_process


def run_model(Model, kin):
    try:
        gen = Model(kin, seed=1)
    except ValueError:
        return None
    return next(gen(1)).final_state()


@pytest.mark.parametrize("frame", ("cms", "ft", "cms2ft", "ft2cms"))
@pytest.mark.parametrize("target", ("p", "air"))
@pytest.mark.parametrize("projectile", ("gamma", "pi-", "p", "He"))
@pytest.mark.parametrize(
    "Model", [m for m in get_all_models() if m not in [im.Sophia20, im.Sibyll21]]
)
def test_final_state_mass_shell(Model, frame, target, projectile):
    p1 = projectile
    p2 = target
    if issubclass(Model, im.EposLHC) and frame == "cms2ft":
        pytest.skip(
            "Epos doesn't conserve mass well when boosted from cms to ft frame."
        )
    if p2 == "air":
        if Model is im.Pythia8:
            pytest.skip("Simulating nuclei in Pythia8 is very time-consuming")

        # cannot use Argon in SIBYLL, so make air from N, O only
        p2 = CompositeTarget((("N", 0.78), ("O", 0.22)))

    if frame == "cms":
        kin = CenterOfMass(100 * GeV, p1, p2)
    elif frame == "ft":
        kin = FixedTarget(TotalEnergy(100 * GeV), p1, p2)
    elif frame == "cms2ft":
        kin = EventKinematicsWithRestframe(
            p1, p2, ecm=100 * GeV, frame=EventFrame.FIXED_TARGET
        )
    elif frame == "ft2cms":
        kin = EventKinematicsWithRestframe(
            p1, p2, elab=100 * GeV, frame=EventFrame.CENTER_OF_MASS
        )
    else:
        assert False  # we should never arrive here

    event = run_in_separate_process(run_model, Model, kin)
    if event is None:
        assert abs(kin.p1) not in Model.projectiles or abs(kin.p2) not in Model.targets
        return
    inv_mass2 = event.en**2 - event.p_tot**2
    inv_mass = np.sign(inv_mass2) * np.sqrt(np.abs(inv_mass2))
    # Increase tolerance above 5 MeV for single-precision models
    atol = 0.005
    if issubclass(Model, im.EposLHC):
        atol = 0.02
    elif Model in [im.QGSJetIII] and p1 == "He" and target == "air":
        atol = 0.05
    assert_allclose(inv_mass, event.m, atol=atol)


@pytest.mark.parametrize("frame", ("cms", "generic"))
@pytest.mark.parametrize(
    "Model", (im.Pythia8, im.Phojet112, im.Phojet191, im.Phojet193)
)
def test_mass_shell_gg(frame, Model):
    p1 = "gamma"
    p2 = "gamma"

    if frame == "cms":
        kin = CenterOfMass(100 * GeV, p1, p2)
    elif frame == "generic":
        kin = EventKinematicsMassless(
            p1, p2, beam=(80, -20), frame=EventFrame.CENTER_OF_MASS
        )
    else:
        assert False  # we should never arrive here

    event = run_in_separate_process(run_model, Model, kin)
    if event is None:
        assert abs(kin.p1) not in Model.projectiles or abs(kin.p2) not in Model.targets
        return
    inv_mass2 = event.en**2 - event.p_tot**2
    inv_mass = np.sign(inv_mass2) * np.sqrt(np.abs(inv_mass2))
    # Increase tolerance above 5 MeV for single-precision models
    atol = 0.005
    assert_allclose(inv_mass, event.m, atol=atol)
