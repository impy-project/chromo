from impy.constants import TeV
from impy.kinematics import EventKinematics
from impy import models
from collections import Counter
import pytest


@pytest.mark.parametrize(
    "model",
    [
        models.Sophia20,
        models.EposLHC,
        models.Sibyll21,
        models.Sibyll23,
        models.Sibyll23c00,
        models.Sibyll23c01,
        models.Sibyll23c02,
        models.Sibyll23c03,
        models.Sibyll23d,
        models.QGSJetII03,
        models.QGSJetII04,
        # models.UrQMD34,
        #   Abort Python with "no collision problem in UrQMD"
        # models.QGSJet01c,
        #   AttributeError: module 'impy.models.qgs01' has no attribute 'qgarr7'
        # models.Pythia6,
        #   Stuck with 100 % CPU
    ],
)
def test_new_interface(model):

    # AF: This is what the user interaction has to yield.
    # It is the typical expected configuration that one
    # wants to run (read pp-mode at energies not exceeding
    # 7 TeV). If you want cosmic ray energies, this should
    # be rather p-N at 10 EeV and lab frame (not yet defined).

    p1pdg = -211  # pi-
    p2pdg = 2212  # proton
    if model is models.Sophia20:
        # Sophia can only do γp, γn
        p1pdg = 22  # gamma
    elif model is models.Pythia6:
        # Pythia6 can only do ee, ep, pp
        p1pdg = 2212  # proton

    ekin = EventKinematics(
        ecm=7 * TeV,
        p1pdg=p1pdg,
        p2pdg=p2pdg,
    )
    gen = model(ekin)

    c = Counter()
    for event in gen(10):
        event.filter_final_state()
        assert len(event.p_ids) > 0
        c.update(event.p_ids)

    assert c[211] > 0, "pi+"
    assert c[-211] > 0, "pi-"
    assert c[2212] > 0, "p"
    assert c[-2212] > 0, "pbar"
