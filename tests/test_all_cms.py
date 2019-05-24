import sys
import os
import numpy as np

root_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(root_dir)
sys.path.append(os.path.join(root_dir, '../DPMJET-III-gitlab'))
sys.path.append(os.path.join(root_dir, '../../apps/pythia8240/lib'))
print root_dir

from impy.definitions import *
from impy.constants import *
from impy.kinematics import EventKinematics
from impy.common import impy_config, pdata

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

impy_config["user_frame"] = 'center-of-mass'

gen_list = ['SIBYLL23C', 'SIBYLL23', 'SIBYLL21', 'DPMJETIII306', 'DPMJETIII171', 
           'EPOSLHC','PHOJET112','PHOJET171','URQMD34','PYTHIA8','QGSJET01C',
           'QGSJETII03','QGSJETII04']

failed = []
passed = []

psrap = {}
eta_range = np.linspace(-5,5,21)
nevents = 5000
for gen in gen_list:
    hist = np.zeros(len(eta_range) - 1) 
    try:
        generator = make_generator_instance(interaction_model_by_tag[gen])
        generator.init_generator(event_kinematics)
        for event in generator.event_generator(event_kinematics, nevents):
            event.filter_final_state_charged()
            hist += 1/float(nevents)*np.histogram(event.eta,bins=eta_range)[0]
        psrap[gen] = hist
    except:
        failed.append(gen)
        continue
        
    passed.append(gen)
    
import IPython
IPython.embed()
