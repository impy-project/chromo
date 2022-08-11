from impy.constants import TeV
from impy.kinematics import EventKinematics
from impy import models
from collections import Counter
import pytest


@pytest.mark.parametrize(
    "model",
    [
        models.EposLHC,
        models.Sibyll21,
        models.Sibyll23d,
        models.QGSJetII03,
        models.QGSJetII04,
        # models.UrQMD34,
        #   Does not compile, see comment in CMakeLists.txt
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

    ekin = EventKinematics(
        ecm=7 * TeV,
        p1pdg=-211,
        # nuc1_prop=(12,6),
        nuc2_prop=(12, 6),
    )

    # # TODO can this be removed?
    # impy_config["user_frame"] = "laboratory"

    gen = model(ekin)

    c = Counter()
    for event in gen(10):
        event.filter_final_state()
        assert len(event.p_ids) > 0
        assert event.impact_parameter < 10
        if model not in (models.Sibyll21, models.Sibyll23d):
            # Sibyll fails these, is this expected?
            assert event.impact_parameter > 0
            assert event.n_wounded_A == 1
        if model not in (
            models.Sibyll21,
            models.Sibyll23d,
            models.QGSJetII03,
            models.QGSJetII04,
        ):
            # Sibyll and QGSJetII fail this, is this expected?
            assert event.n_wounded_B > 0
        # assert event.n_NN_interactions > 0

        c.update(event.p_ids)

    assert c[211] > 0, "pi+"
    assert c[-211] > 0, "pi-"
    assert c[2212] > 0, "p"
    assert c[-2212] > 0, "pbar"
