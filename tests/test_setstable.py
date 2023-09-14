import numpy as np
from chromo.constants import GeV
from chromo.kinematics import CenterOfMass
from .util import run_in_separate_process
from chromo.util import get_all_models
import pytest
from collections import Counter
from particle import literals as lp
import platform


decay_list = [
    lp.pi_plus.pdgid,
    lp.K_plus.pdgid,
    lp.pi_0.pdgid,
    lp.K_S_0.pdgid,
    lp.K_L_0.pdgid,
    lp.Lambda.pdgid,
]


def run_model(Model, stable):
    p1 = "p"
    if Model.name == "Sophia":
        p1 = "gamma"
    p2 = "p"
    kin = CenterOfMass(200 * GeV, p1, p2)
    model = Model(kin, seed=1)
    for pid in decay_list:
        model.set_stable(pid, stable)
    c = Counter()
    for event in model(100):
        ev = event.final_state()
        c.update(ev.pid)
    return c


@pytest.mark.parametrize("stable", (False, True))
@pytest.mark.parametrize("Model", get_all_models())
def test_setstable(Model, stable):
    if "UrQMD" in Model.name:
        pytest.xfail(f"{Model.pyname} does not support changing decays")
    if not stable:
        if ("QGSJet" in Model.name) and (platform.system() == "Windows"):
            pytest.xfail(f"{Model.pyname} does not support changing decays")

    c = run_in_separate_process(run_model, Model, stable)

    csel = np.array([c[pid] for pid in decay_list])
    if stable:
        assert np.all(csel > 0), f"stable={stable} counters must be > 0, got {csel}"
    else:
        assert np.all(csel == 0), f"stable={stable} counters must be 0, got {csel}"
