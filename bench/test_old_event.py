from impy import models as im
from impy.constants import TeV
from impy.kinematics import EventKinematics
import pytest

N = 1000
ekin = EventKinematics(ecm=10 * TeV, particle1=2212, particle2=2212)
models = {name: getattr(im, name)(ekin) for name in ("Sibyll21", "Pythia6")}


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
            event.filter_final_state()

    benchmark(run, model)


@pytest.mark.parametrize("model", models)
def test_event_final_charged(model, benchmark):
    model = models[model]

    def run(model):
        for event in model(N):
            event.filter_final_state_charged()

    benchmark(run, model)
