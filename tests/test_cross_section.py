import pytest
from impy.util import get_all_models
from impy.kinematics import CompositeTarget, CenterOfMass, TeV
from impy import models as im


@pytest.mark.parametrize("Model", get_all_models())
@pytest.mark.parametrize("target", ("p", "O", "air"))
@pytest.mark.parametrize("projectile", ("gamma", "pi-", "p", "He"))
def test_cross_section(Model, projectile, target):
    p1 = projectile
    p2 = target
    if p2 == "air":
        if Model is im.Pythia8:
            pytest.skip("Simulating nuclei in Pythia8 is very time-consuming")

        # cannot use Argon in SIBYLL, so make air from N, O only
        p2 = CompositeTarget((("N", 0.8), ("O", 0.2)))

    kin = CenterOfMass(1 * TeV, p1, p2)

    model = Model(seed=1)

    if abs(kin.p1) not in Model.projectiles or abs(kin.p2) not in Model.targets:
        with pytest.raises(ValueError):
            model.cross_section(kin)
        return

    # known issues
    pyname = Model.pyname
    if pyname.startswith("Sibyll") and (kin.p2.A or 1) > 1:
        pytest.xfail("h-A should work, but need to be fixed in Sibyll")
    if pyname == "QGSJet01d" and (kin.p1.A or 1) > 1:
        pytest.xfail("no support for nuclear projectiles in QGSJet01d")
    if pyname == "Pythia8" and ((kin.p1.A or 1) > 1 or (kin.p2.A or 1) > 1):
        pytest.xfail("A-A and h-A cross-sections are not yet supported for Pythia8")
    if pyname == "DpmjetIII191" and projectile == "He" and target == "air":
        pytest.xfail("DpmjetIII191 returns NaN for (He, air)")
    if any(
        [
            (pyname == "DpmjetIII193" and kin.p1.is_nucleus and kin.p2.is_nucleus),
            (pyname.startswith("Dpmjet") and abs(kin.p1) == 211),
        ]
    ):
        pytest.xfail(f"{Model.pyname} aborts on ({projectile}, {target})")

    c = model.cross_section(kin)

    if Model is im.UrQMD34:
        # only provides total cross-section
        # HD: these values are "accepted, but make no sense
        expected = {
            ("pi-", "p"): None,
            ("pi-", "air"): 1200,
            ("p", "p"): None,
            ("p", "air"): 1200,
            ("He", "p"): 900,
            ("He", "air"): 1400,
            ("pi-", "O"): 1266,
            ("He", "O"): 1500,
            ("p", "O"): 1200,
        }[(projectile, target)]

        if expected is None:
            pytest.xfail("UrQMD34 returns 0 for p p and pi- p")

        assert c.total == pytest.approx(expected, rel=0.1)
        return

    # These are "average" values derived from the output
    # of several models. They just indicate that the
    # model gives an answer in the expected ballpark.
    #
    # Models seem to agree on p-p, so we base other
    # reference values on that.
    # - p-A should be roughly A^(2/3) * p-p
    # - pi-p should be roughly (2/3)^(2/3) * p-p
    # - gamma-p should be 1/137 * (1/3)^(2/3), assuming alpha=1
    pp = 53
    power = 2 / 3  # volume to area
    f_pi = (2 / 3) ** power
    f_O = 16**power
    f_He = 4**power
    f_air = 0.8 * 14**power + 0.2 * f_O
    f_gamma = (1 / 3) ** power / 137
    expected = {
        ("gamma", "p"): pp * f_gamma,
        ("pi-", "p"): pp * f_pi,
        ("pi-", "O"): pp * f_pi * f_O,
        ("pi-", "air"): pp * f_pi * f_air,
        ("p", "p"): pp,
        ("p", "O"): pp * f_O,
        ("p", "air"): pp * f_air,
        ("He", "p"): pp * f_He,
        ("He", "air"): pp * f_He * f_air,
        ("He", "O"): pp * f_He * f_O,
    }[(projectile, target)]

    # large tolerance to accommodate substantial differences between models
    assert c.inelastic == pytest.approx(expected, rel=0.6)
