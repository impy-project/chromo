# FIXME this is redundant with test_generators.py
# instead do specific tests that only work for eposlhc here

from impy.constants import GeV
from impy.kinematics import EventKinematics
from impy import impy_config
from impy.models import EposLHC
from collections import Counter


def test_eposlhc():
    event_kinematics = EventKinematics(
        ecm=7000 * GeV,
        particle1=2212,
        particle2=2212
        # particle2=(14,7)
    )

    # FIXME is this still needed
    impy_config["user_frame"] = "laboratory"

    generator = EposLHC(event_kinematics)

    c = Counter()
    for event in generator(2):
        # generator.lib.pydat3.mdcy[102 - 1, 0] = 1
        ev = event.final_state_charged()
        c.update(ev.pid)

    assert c[211] > 0
    assert c[2212] > 0
