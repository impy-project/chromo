import sys
import os
sys.path.append(os.path.dirname(__file__) + "/..")

import numpy as np

from impy.constants import *
from impy.kinematics import EventKinematics
from impy.models.sibyll import SIBYLLRun


eventkinematics = EventKinematics(ecm=7 * TeV, p1pdg=2212, p2pdg=2122)


# Something like this can be a template for passing args
config = dict(
    label='SIBYLL2.3c_testrun',
    event_config={'charged_only': True},
    event_kinematics=eventkinematics,
)

# This is a workaround. Something needs to keep track of the
# fortran imports.
import sib23c

generator = SIBYLLRun(sib23c, **config)

generator.init_generator()

for event in generator.event_generator(eventkinematics, 50):
    print event.px, event.py, event.pz, event.en