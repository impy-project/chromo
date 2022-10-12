# FIXME this is redundant with test_generators.py
# instead do specific tests that only work for eposlhc here

from impy.constants import TeV
from impy.kinematics import EventKinematics
from impy.models import Phojet112
from collections import Counter


def test_phojet():

    # AF: This is what the user interaction has to yield.
    # It is the typical expected configuration that one
    # wants to run (read pp-mode at energies not exceeding
    # 7 TeV). If you want cosmic ray energies, this should
    # be rather p-N at 10 EeV and lab frame (not yet defined).

    p1pdg = -211  # pi-
    p2pdg = 2212  # proton

    evt_kin = EventKinematics(
        ecm=7 * TeV,
        particle1=p1pdg,
        particle2=p2pdg,
    )

    # Some models need to initialize same fortran code,
    # which can only be initialized once, therefore run
    # in separate thread

    gen = Phojet112(evt_kin)

    c = Counter()
    for event in gen(10):
        event.filter_final_state()
        assert len(event.pid) > 0
        c.update(event.pid)

    assert c[211] > 0, "pi+"
    assert c[-211] > 0, "pi-"
    assert c[2212] > 0, "p"
    assert c[-2212] > 0, "pbar"
