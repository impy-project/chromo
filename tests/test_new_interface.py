
import sys
import os
import numpy as np

sys.path.append(os.path.dirname(__file__) + "/..")

from impy.constants import *
from impy.kinematics import EventKinematics
from impy.models.sibyll import SIBYLLRun
from impy.models.dpmjetIII import DpmjetIIIRun
from impy.common import impy_config, pdata




# AF: This is what the user interaction has to yield.
# It is the typical expected configuration that one
# wants to run (read pp-mode at energies not exceeding
# 7 TeV). If you want cosmic ray energies, this should
# be rather p-N at 10 EeV and lab frame (not yet defined).

event_kinematics = EventKinematics(
    ecm=7 * TeV, 
    p1pdg=2212, 
    # nuc1_prop=(12,6),
    nuc2_prop=(16,8))

impy_config["user_frame"] = 'laboratory'

libhandle = None
# Run legacy DPMJET
exec 'import {0} as libhandle'.format('dpmjet306')
generator = DpmjetIIIRun(libhandle) 
# exec 'import {0} as libhandle'.format('sib23c')
# generator = SIBYLLRun(libhandle)

# If init remains without args, it should go to the contructor.
generator.init_generator(event_kinematics)

# This  
for event in generator.event_generator(event_kinematics, 50):
    # print 'px', event.px
    # print 'py', event.py 
    print 'pz', event.pz
    print 'en', event.en
    print 'p_ids', event.p_ids
    # print event.impact_parameter, event.n_wounded_A, event.n_wounded_B#, event.n_NN_interactions 

# AF: Maybe it would be better, or a good alternative
# to make a Particle class and the MCEvent provides an iterator
# to it's contents as objects of this Particle class. This would
# be Pythia 8 style but probably slow in PYTHON.

# for event in ....:
#     for particle in event:
#         print particle.p_id, particle.en...

