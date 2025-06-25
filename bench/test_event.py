import pytest

from chromo import models as im
from chromo.constants import TeV
from chromo.kinematics import CenterOfMass

N = 1000
evt_kin = CenterOfMass(10 * TeV, "p", "p")
models = {name: getattr(im, name)(evt_kin) for name in ("Sibyll21", "Pythia6")}


@pytest.mark.parametrize("model", models)
def test_event(model, benchmark):
    model = models[model]

    def run(model):
        for event in model(N):
            pass

    benchmark(run, model)


@pytest.mark.parametrize("model", models)
def test_event_final(model, benchmark):
    model = models[model]

    def run(model):
        for event in model(N):
            event.final_state()

    benchmark(run, model)


@pytest.mark.parametrize("model", models)
def test_event_final_charged(model, benchmark):
    model = models[model]

    def run(model):
        for event in model(N):
            event.final_state_charged()

    benchmark(run, model)
