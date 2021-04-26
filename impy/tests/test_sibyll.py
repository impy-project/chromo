import sys
import os
import numpy as np

root_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(root_dir)
sys.path.append(os.path.join(root_dir, '../DPMJET-III-gitlab'))

from impy.definitions import *
from impy.constants import *
from impy.kinematics import EventKinematics
from impy import impy_config, pdata

# AF: This is what the user interaction has to yield.
# It is the typical expected configuration that one
# wants to run (read pp-mode at energies not exceeding
# 7 TeV). If you want cosmic ray energies, this should
# be rather p-N at 10 EeV and lab frame (not yet defined).

event_kinematics = EventKinematics(ecm=7000 * GeV,
                                   p1pdg=2212,
                                   p2pdg=2212
                                   # nuc2_prop=(14,7)
                                   )

impy_config["user_frame"] = 'laboratory'

generator = make_generator_instance(interaction_model_by_tag['SIBYLL23C'])
generator.init_generator(event_kinematics)
# import IPython
# IPython.embed()

# This
for event in generator.event_generator(event_kinematics, 2):
    # generator.lib.pydat3.mdcy[102 - 1, 0] = 1
    import IPython
    IPython.embed()
    event.filter_final_state_charged()
    # print 'px', event.px
    # print 'py', event.py
    # print 'pz', event.pz
    # print 'en', event.en
    print 'p_ids', event.p_ids
    # print 'impact param', event.impact_parameter
