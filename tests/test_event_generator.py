from impy.models import EposLHC, Sibyll21
from impy.kinematics import EventKinematics
from impy.constants import TeV
import numpy as np
import pytest

proton_pid = 2212


@pytest.mark.parametrize("model", (EposLHC, Sibyll21))
def test_model(model):
    ekin = EventKinematics(ecm=13 * TeV, p1pdg=proton_pid, p2pdg=proton_pid)
    gen = model(ekin)
    for event in gen(100):
        assert np.sum(event.px)
        assert np.sum(event.py)
        assert np.sum(event.pz)
        assert np.sum(event.en)
        assert np.sum(event.p_ids)
