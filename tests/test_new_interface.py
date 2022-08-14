from impy.constants import TeV
from impy.kinematics import EventKinematics
import impy.models as im
from impy.models import Sophia20
from impy.common import MCRun
from collections import Counter
import pytest
from multiprocessing import Pool


models = set()
for key in dir(im):
    obj = getattr(im, key)
    if hasattr(obj, "__mro__") and MCRun in obj.__mro__:
        models.add(obj)


def run_model(model, ekin):
    gen = model(ekin)

    c = Counter()
    for event in gen(10):
        event.filter_final_state()
        assert len(event.p_ids) > 0
        c.update(event.p_ids)

    return c


@pytest.mark.parametrize("model", models)
def test_new_interface(model):

    # AF: This is what the user interaction has to yield.
    # It is the typical expected configuration that one
    # wants to run (read pp-mode at energies not exceeding
    # 7 TeV). If you want cosmic ray energies, this should
    # be rather p-N at 10 EeV and lab frame (not yet defined).

    p1pdg = -211  # pi-
    p2pdg = 2212  # proton
    if model is Sophia20:
        # Sophia can only do γp, γn
        p1pdg = 22  # gamma
    # elif model is models.Pythia6:
    #     # Pythia6 can only do ee, ep, pp
    #     p1pdg = 2212  # proton

    ekin = EventKinematics(
        ecm=7 * TeV,
        p1pdg=p1pdg,
        p2pdg=p2pdg,
    )

    # Some models need to initialize same fortran code,
    # which can only be initialized once, therefore run
    # in separate thread
    with Pool(1) as p:
        r = p.apply_async(run_model, (model, ekin))
        try:
            c = r.get(timeout=15)
        except TimeoutError:
            assert False

    assert c[211] > 0, "pi+"
    assert c[-211] > 0, "pi-"
    assert c[2212] > 0, "p"
    assert c[-2212] > 0, "pbar"
