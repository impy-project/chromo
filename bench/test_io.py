import os
from impy import models as im
from impy.constants import TeV
from impy.kinematics import EventKinematics
import pytest
import pyhepmc

N = 1000
ekin = EventKinematics(ecm=10 * TeV, p1pdg=2212, p2pdg=2212)
models = {name: getattr(im, name)(ekin) for name in ("Sibyll21", "Pythia6")}


@pytest.mark.parametrize("model", models)
def test_event(model, benchmark):
    model = models[model]

    def run(model):
        with pyhepmc.open("output.dat", "w") as f:
            for event in model(N):
                f.write(event)

    benchmark(run, model)

    os.unlink("output.dat")
