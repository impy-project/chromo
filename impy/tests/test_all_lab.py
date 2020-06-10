from __future__ import print_function

import sys
import os
import numpy as np
from multiprocessing import Pool
import tempfile
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

event_kinematics = EventKinematics(ecm=7000 * GeV,
                                   p1pdg=2212,
                                   p2pdg=2212
                                   # nuc2_prop=(14,7)
                                   )

impy_config["user_frame"] = 'center-of-mass'

gen_list = [
    'SIBYLL23C', 
    'SIBYLL23', 
    'SIBYLL21', 
    'DPMJETIII306', 
    'DPMJETIII171', 
    'EPOSLHC',
#     'PHOJET112',
#     'PHOJET171', # Does not support nuclei
    'URQMD34',
    'PYTHIA8',
    'QGSJET01C',
    'QGSJETII03',
    'QGSJETII04'
]

xlab_bins = np.linspace(0,1,21)
xlab_widths = xlab_bins[1:] - xlab_bins[:-1]
xlab_centers = 0.5*(xlab_bins[1:] + xlab_bins[:-1])
nevents = 5000
norm = 1./float(nevents)/xlab_widths

def run_generator(gen,*args):
    print('Testing',gen)
    hist_p = np.zeros(len(xlab_centers))
    hist_pi = np.zeros(len(xlab_centers)) 
    try:
        log = tempfile.mkstemp()[1]
        generator = make_generator_instance(interaction_model_by_tag[gen])
        generator.init_generator(event_kinematics,logfname=log)
        for event in generator.event_generator(event_kinematics, nevents):
            event.filter_final_state_charged()

            hist_p += np.histogram(event.xlab[event.p_ids == 2212],
                               bins=xlab_bins,
                               weights=event.xlab[event.p_ids == 2212]**1.7)[0]

            hist_pi += np.histogram(event.xlab[np.abs(event.p_ids) == 211],
                                bins=xlab_bins,
                                weights=event.xlab[np.abs(event.p_ids) == 211]**1.7)[0]
            
        return True, gen, log, hist_p, hist_pi
    except:
        return False, gen, log, hist_p, hist_pi
        
pool = Pool(processes=8)
result = [pool.apply_async(run_generator, (gen,)) for gen in gen_list]
result = [res.get(timeout=100000) for res in result]

logs = {}
xlab_protons = {}
xlab_piplus = {}
failed = []
passed = []

for r, gen, log, hist_p, hist_pi in result:
    if r:
        passed.append(gen)
        xlab_protons[gen] = hist_p
        xlab_piplus[gen] = hist_pi
    else:
        failed.append(gen)
        
    with open(log) as f:
            logs[gen] = f.read()
        

info(0, 'Test results for 158 GeV pC collisions in lab frame:\n')
info(0, 'Passed:', '\n', '\n '.join(passed))
info(0, '\nFailed:', '\n', '\n '.join(failed))

import pickle
pickle.dump((xlab_bins, xlab_protons, xlab_piplus, logs),
            open(os.path.splitext(__file__)[0] + '.pkl','wb'), protocol=-1)
