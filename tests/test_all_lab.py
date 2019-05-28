from __future__ import print_function

import sys
import os
import numpy as np

root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))
sys.path.append(root_dir)
sys.path.append(os.path.join(root_dir, '../DPMJET-III-gitlab'))
sys.path.append(os.path.join(root_dir, '../../apps/pythia8240/lib'))

from impy.definitions import *
from impy.constants import *
from impy.kinematics import EventKinematics
from impy.common import impy_config, pdata
from impy.util import info

# AF: This is what the user interaction has to yield.
# It is the typical expected configuration that one
# wants to run (read pp-mode at energies not exceeding
# 7 TeV). If you want cosmic ray energies, this should
# be rather p-N at 10 EeV and lab frame (not yet defined).

event_kinematics = EventKinematics(plab=158 * GeV,
                                   p1pdg=2212,
                                   nuc2_prop=(12,6)
                                   )

impy_config["user_frame"] = 'laboratory'

gen_list = ['SIBYLL23C', 'SIBYLL23', 'SIBYLL21', 'DPMJETIII306', 'DPMJETIII171', 
           'EPOSLHC','PHOJET112','PHOJET171','URQMD34','PYTHIA8','QGSJET01C',
           'QGSJETII03','QGSJETII04']

failed = []
passed = []

xlab_hist = {}
xlab_range = np.linspace(0,1,21)
nevents = 5000
for gen in gen_list:
    hist = np.zeros(len(xlab_range) - 1) 
    try:
        generator = make_generator_instance(interaction_model_by_tag[gen])
        generator.init_generator(event_kinematics)
        for event in generator.event_generator(event_kinematics, nevents):
            event.filter_final_state_charged()
            hist += 1/float(nevents)*np.histogram(event.xlab,bins=xlab_range,
                                                  weights=event.xlab**1.7)[0]
        xlab_hist[gen] = hist
    except:
        failed.append(gen)
        continue
        
    passed.append(gen)
  
info(0, 'Test results for 158 GeV pC collisions in lab frame:\n')
info(0, 'Passed:', '\n', '\n '.join(passed))
info(0, '\nFailed:', '\n', '\n '.join(failed))

import pickle
pickle.dump((xlab_range, xlab_hist),
            open(os.path.splitext(__file__)[0] + '.pkl','wb'), protocol=-1)

# import IPython
# IPython.embed()
