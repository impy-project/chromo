import numpy as np

from impy.models import QGSJet01c, EposLHC
from impy.constants import GeV
from impy.kinematics import EventKinematics
from impy import impy_config
import pytest
from collections import Counter

# TODO this test is broken


@pytest.mark.parametrize("model", [QGSJet01c, EposLHC])
def test_unstable(model):
    event_kinematics = EventKinematics(ecm=200 * GeV, particle1=2212, particle2=2212)

    # impy_config["user_frame"] = 'laboratory'
    impy_config["tau_stable"] = np.inf
    # impy_config['pi0_stable'] = False
    generator = model(event_kinematics)

    # decay_list = [211, 321, 111, 2112, 310, 130, 3122, -3212]
    decay_list = [211, 321, 111, 2112, 310, 130]

    for stable in [True, False]:
        for pid in decay_list:
            generator.set_stable(pid, stable)

        c = Counter()
        for event in generator(100):
            # generator.lib.pydat3.mdcy[102 - 1, 0] = 1
            ev = event.final_state()

            c.update(ev.pid)

            # for pid in decay_list:
            #     if pid in np.abs(event.pid):
            #         print('Decay not working for',pid)
            # raise Exception('Decay not working for',pid)
            # print(event.pid)
            # print 'px', event.px
            # print 'py', event.py
            # print 'pz', event.pz
            # print 'en', event.en
            # print 'pid', event.pid
            # print 'impact param', event.impact_parameter

        if stable:
            for pid in decay_list:
                assert c[pid] > 0, f"stable={stable} pid={pid}: {c[pid]} > 0 violated"
        else:
            for pid in decay_list:
                assert c[pid] == 0, f"stable={stable} pid={pid}: {c[pid]} == 0 violated"
