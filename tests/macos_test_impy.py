import impy
import impy.models as im
import pytest


models = [im.Sophia20, im.Sibyll21, im.Sibyll23, im.Sibyll23c, im.Sibyll23d]


@pytest.mark.parametrize("model", models)
def test_generator(model):
    print(model.label)
    if model == im.Sophia20:
        p1 = "photon"
    else:
        p1 = "proton"

    ekin = impy.kinematics.FixedTarget(1000, p1, "proton")
    gen = model(ekin)

    for event in gen(1):
        assert len(event.pid) > 0
