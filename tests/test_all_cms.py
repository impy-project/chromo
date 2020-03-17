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
from impy import impy_config, pdata
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
    'DPMJETIII191', 
    'EPOSLHC',
    'PHOJET112',
    'PHOJET191',
    'URQMD34',
    'PYTHIA8',
    'QGSJET01C',
    'QGSJETII03',
    'QGSJETII04'
]

eta_bins = np.linspace(-5,5,21)
eta_widths = eta_bins[1:] - eta_bins[:-1]
eta_centers = 0.5*(eta_bins[1:] + eta_bins[:-1])
nevents = 5000
norm = 1./float(nevents)/eta_widths

def run_generator(gen,*args):
    print('Testing',gen)
    hist = np.zeros(len(eta_centers)) 
    try:
        log = tempfile.mkstemp()[1]
        generator = make_generator_instance(interaction_model_by_tag[gen])
        generator.init_generator(event_kinematics,logfname=log)
        for event in generator.event_generator(event_kinematics, nevents):
            event.filter_final_state_charged()
            hist += norm*np.histogram(event.eta,bins=eta_bins)[0]
        return True, gen, log, hist
    except:
        return False, gen, log, hist
        
pool = Pool(processes=8)
result = [pool.apply_async(run_generator, (gen,)) for gen in gen_list]
result = [res.get(timeout=100000) for res in result]

failed = []
passed = []
psrap = {}
logs = {}

for r, gen, log, hist in result:
    if r:
        passed.append(gen)
        psrap[gen] = hist
    else:
        failed.append(gen)
        
    with open(log) as f:
            logs[gen] = f.read()
        
info(0, 'Test results for 7 TeV pp collisions in cms frame:\n')
info(0, 'Passed:', '\n', '\n '.join(passed))
info(0, '\nFailed:', '\n', '\n '.join(failed))

import pickle
pickle.dump((eta_bins, psrap, logs),
            open(os.path.splitext(__file__)[0] + '.pkl','wb'), protocol=-1)

# import IPython
# IPython.embed()
