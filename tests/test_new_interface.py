import sys
import os
import numpy as np

root_dir = os.path.abspath(os.path.dirname(__file__) + "/..")
sys.path.append(root_dir)
sys.path.append(os.path.join(root_dir,'../DPMJET-III-gitlab'))


from impy.definitions import *
from impy.constants import *
from impy.kinematics import EventKinematics
from impy.common import impy_config, pdata

# AF: This is what the user interaction has to yield.
# It is the typical expected configuration that one
# wants to run (read pp-mode at energies not exceeding
# 7 TeV). If you want cosmic ray energies, this should
# be rather p-N at 10 EeV and lab frame (not yet defined).

event_kinematics = EventKinematics(
    ecm=7 * TeV,
    p1pdg=-211,
    # nuc1_prop=(12,6),
    nuc2_prop=(12, 6))

impy_config["user_frame"] = 'laboratory'

generator = make_generator_instance(interaction_model_by_tag['DPMJETIII171'])
generator.init_generator(event_kinematics)
# import IPython
# IPython.embed()

# This
for event in generator.event_generator(event_kinematics, 10):
    event.filter_final_state()
    # print 'px', event.px
    # print 'py', event.py
    # print 'pz', event.pz
    # print 'en', event.en
    print('p_ids', event.p_id)
    print('impact param', event.impact_parameter)
    # import IPython
    # IPython.embed()
    # print event.impact_parameter, event.n_wounded_A, event.n_wounded_B#, event.n_NN_interactions
