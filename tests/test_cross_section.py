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
        p2 = CompositeTarget((("N", 0.78), ("O", 0.22)))

    kin = CenterOfMass(1 * TeV, p1, p2)

    model = Model(seed=1)

    if abs(kin.p1) not in Model.projectiles or abs(kin.p2) not in Model.targets:
        with pytest.raises(ValueError):
            model.cross_section(kin)
        return

    # known issues
    if Model.name == "DPMJET-III" and not kin.p1.is_nucleus and kin.p2.is_nucleus:
        pytest.xfail("DPMJet rejects h-A, although it can do h-p, p-A, and A-A")
    if Model.name == "SIBYLL" and (kin.p2.A or 1) > 1:
        pytest.xfail("h-A should work, but need to be fixed in Sibyll")
    if Model.pyname == "QGSJet01d" and (kin.p1.A or 1) > 1:
        pytest.xfail("no support for nuclear projectiles in QGSJet01d")
    if Model.pyname == "Pythia8" and (kin.p1.A or 1) > 1 or (kin.p2.A or 1) > 1:
        pytest.xfail("A-A and h-A cross-sections are not yet supported for Pythia8")

    c = model.cross_section(kin)

    if Model is im.UrQMD34:
        # only provides total cross-section
        expected = {
            ("pi-", "p"): 0,
            ("pi-", "air"): 1200,
            ("p", "p"): 0,
            ("p", "air"): 1200,
            ("He", "p"): 900,
            ("He", "air"): 1400,
        }[(projectile, target)]

        assert c.total == pytest.approx(expected, rel=0.1)
        return

    # p-A should be roughly A^(2/3) * p-p
    # pi-p should be roughly (2/3)^(2/3) * p-p
    expected = {
        ("gamma", "p"): 0.21,
        ("pi-", "p"): 40,
        ("pi-", "O"): 180,
        ("pi-", "air"): 150,
        ("p", "p"): 53,
        ("p", "O"): 240,
        ("p", "air"): 220,
        ("He", "p"): 130,
        ("He", "air"): 350,
        ("He", "O"): 380,
    }[(projectile, target)]

    assert c.inelastic == pytest.approx(expected, rel=0.3)
