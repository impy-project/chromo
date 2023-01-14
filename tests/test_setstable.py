import numpy as np
from impy.constants import GeV
from impy.kinematics import CenterOfMass
from impy.util import get_all_models
import pytest
from collections import Counter
from particle import literals as lp


DECAY_LIST = [
    lp.pi_plus.pdgid,
    lp.K_plus.pdgid,
    lp.pi_0.pdgid,
    lp.K_S_0.pdgid,
    lp.K_L_0.pdgid,
]


def run_model(Model, stable):
    p1 = "p"
    if Model.name == "Sophia":
        p1 = "gamma"
    p2 = "p"

    model = Model(seed=1)
    for pid in DECAY_LIST:
        model.set_stable(pid, stable)

    c = Counter()
    kin = CenterOfMass(200 * GeV, p1, p2)
    for event in model(kin, 100):
        ev = event.final_state()
        c.update(ev.pid)
    return c


@pytest.mark.parametrize("stable", (False, True))
@pytest.mark.parametrize("Model", get_all_models())
def test_setstable(Model, stable):
    if not stable:
        if any(part in Model.name for part in ("DPMJET", "PhoJet")):
            pytest.xfail(
                f"{Model.pyname} does not support decaying charged pions or any kaons"
            )
        if "QGSJet" in Model.name:
            pytest.xfail(f"{Model.pyname} does not support changing decays")
    if "UrQMD" in Model.name:
        pytest.xfail(f"{Model.pyname} does not support changing decays")

    c = run_model(Model, stable)

    csel = np.array([c[pid] for pid in DECAY_LIST])
    if stable:
        assert np.all(csel > 0), f"stable={stable} counters must be > 0, got {csel}"
    else:
        assert np.all(csel == 0), f"stable={stable} counters must be 0, got {csel}"
