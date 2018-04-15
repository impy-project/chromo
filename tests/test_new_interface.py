
import sys
import os
import numpy as np

from impy.constants import *
from impy.kinematics import EventKinematics
from impy.models.sibyll import SIBYLLRun

sys.path.append(os.path.dirname(__file__) + "/..")


# AF: This is what the user interaction has to yield.
# It is the typical expected configuration that one
# wants to run (read pp-mode at energies not exceeding
# 7 TeV). If you want cosmic ray energies, this should
# be rather p-N at 10 EeV and lab frame (not yet defined).

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

# Create some wrapper around the fortran library
generator = SIBYLLRun(sib23c, **config)
# If init remains without args, it should go to the contructor.
generator.init_generator()

# This 
for event in generator.event_generator(eventkinematics, 50):
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

