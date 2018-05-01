
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

event_kinematics = EventKinematics(ecm=7 * TeV, p1pdg=2212, p2pdg=2122)

# Something like this can be a template for passing args
# config = dict(
#     label='SIBYLL2.3c_testrun',
#     event_config={'charged_only': True},
#     event_kinematics=eventkinematics,
# )

# This statement can be replaced by a dynamical exec statement
libhandle = None
exec 'import {0} as libhandle'.format('dpmjet306')
# exec 'import {0} as libhandle'.format('sib23c')
# Create some wrapper around the fortran library
generator = DpmjetIIIRun(libhandle) 
# generator = SIBYLLRun(libhandle, **config)
# If init remains without args, it should go to the contructor.
generator.init_generator(event_kinematics)

# This  
for event in generator.event_generator(event_kinematics, 50):
    print 'px', event.px
    print 'py', event.py 
    print 'pz', event.pz
    print 'en', event.en

# AF: Maybe it would be better, or a good, an alternative
# to make a Particle class and the MCEvent provides an iterator
# to it's contents as objects of this Particle class. This would
# be Pythia 8 style.

# for event in ....:
#     for particle in event:
#         print particle.p_id, particle.en...

