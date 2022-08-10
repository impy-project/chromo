from impy.constants import TeV
from impy.kinematics import EventKinematics
from impy.common import impy_config
from impy.models import EposLHC, DpmjetIII306, QGSJet01c, QGSJetII04  # noqa
from collections import Counter
import pytest


@pytest.mark.parametrize(
    "model",
    [
        EposLHC,
        DpmjetIII306,
        # QGSJet01c, # broken: AttributeError: module 'qgs01' has no attribute 'qgarr7'
        # QGSJetII04, # broken because apparently cannot find some files
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

    # TODO can this be removed?
    impy_config["user_frame"] = "laboratory"

    gen = model(ekin)

    c = Counter()
    for event in gen(10):
        event.filter_final_state()
        # print 'px', event.px
        # print 'py', event.py
        # print 'pz', event.pz
        # print 'en', event.en
        assert len(event.p_ids) > 0
        assert event.impact_parameter > 0
        assert event.impact_parameter < 10
        assert event.n_wounded_A == 1
        assert event.n_wounded_B > 0
        # assert event.n_NN_interactions > 0

        c.update(event.p_ids)

    assert c[211] > 0, "pi+"
    assert c[-211] > 0, "pi-"
    assert c[2212] > 0, "p"
    assert c[-2212] > 0, "pbar"
