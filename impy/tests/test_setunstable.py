import sys
import os
import numpy as np

from impy.definitions import *
from impy.constants import *
from impy.kinematics import EventKinematics
from impy import impy_config, pdata


def test_unstable(gen_name='SIBYLL23D'):
    event_kinematics = EventKinematics(
        ecm=200 * GeV,
        p1pdg=2212,
        p2pdg=2212)

    # impy_config["user_frame"] = 'laboratory'
    impy_config['tau_stable'] = np.inf
    # impy_config['pi0_stable'] = False
    generator = make_generator_instance(interaction_model_by_tag[gen_name])
    generator.init_generator(event_kinematics)

    decay_list = [211,321,111,2112, 310,130,13,-13,3122,-3212]
    for pid in decay_list:
        generator.set_stable(pid, stable=False)
    print(decay_list)
    for ie, event in enumerate(generator.event_generator(event_kinematics, 1000)):
        # generator.lib.pydat3.mdcy[102 - 1, 0] = 1
        event.filter_final_state()
        print(ie)
        # for pid in decay_list:
        #     if pid in np.abs(event.p_ids):
        #         print('Decay not working for',pid)
                # raise Exception('Decay not working for',pid)
        # print(event.p_ids)
        # print 'px', event.px
        # print 'py', event.py
        # print 'pz', event.pz
        # print 'en', event.en
        # print 'p_ids', event.p_ids
        # print 'impact param', event.impact_parameter

if __name__ in ['__main__', '__test__']:
    test_unstable('DPMJETIII191')