from __future__ import print_function

import sys
import os
import numpy as np
from multiprocessing import Pool, freeze_support
import tempfile

root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(root_dir, "../../apps/pythia8240/lib"))

from impy.constants import GeV
from impy.kinematics import EventKinematics
from impy import impy_config
from impy.util import info
from impy.models import (
    Sibyll23d,
    Sibyll23c,
    Sibyll23c01,
    Sibyll23c00,
    Sibyll23,
    Sibyll21,
    DpmjetIII306,
    DpmjetIII191,
    EposLHC,
    Phojet112,
    Phojet191,
    UrQMD34,
    QGSJet01c,
    QGSJetII03,
    QGSJetII04,
)


def run_generator(model):
    event_kinematics = EventKinematics(
        ecm=7000 * GeV,
        p1pdg=2212,
        p2pdg=2212
        # nuc2_prop=(14,7)
    )

    impy_config["user_frame"] = "center-of-mass"

    eta_bins = np.linspace(-5, 5, 21)
    eta_widths = eta_bins[1:] - eta_bins[:-1]
    eta_centers = 0.5 * (eta_bins[1:] + eta_bins[:-1])
    nevents = 5000
    norm = 1.0 / float(nevents) / eta_widths
    hist = np.zeros(len(eta_centers))
    try:
        log = tempfile.mkstemp()[1]
        generator = model(event_kinematics)
        for event in generator(nevents):
            event.filter_final_state_charged()
            hist += np.histogram(event.eta, bins=eta_bins)[0]
        return True, generator.__class__.__name__, log, eta_bins, hist * norm
    except Exception:
        return False, generator.__class__.__name__, log, eta_bins, hist * norm


models = [
    Sibyll23d,
    Sibyll23c,
    Sibyll23c01,
    Sibyll23c00,
    Sibyll23,
    Sibyll21,
    DpmjetIII306,
    DpmjetIII191,
    EposLHC,
    Phojet112,
    Phojet191,
    UrQMD34,
    # 'PYTHIA8',
    QGSJet01c,
    QGSJetII03,
    QGSJetII04,
]


def test_all_cms():
    freeze_support()
    pool = Pool(processes=32)
    result = [pool.apply_async(run_generator, (model,)) for model in models]
    result = [res.get(timeout=1000) for res in result]

    failed = set()
    passed = set()
    psrap = {}
    logs = {}

    eta_bins = result[0][3]

    for r, gen, log, _, hist in result:
        if r:
            passed.add(gen)
            psrap[gen] = hist
        else:
            failed.add(gen)

        with open(log) as f:
            logs[gen] = f.read()

    info(0, "Test results for 7 TeV pp collisions in cms frame:\n")
    info(0, "Passed:", "\n", "\n ".join(passed))
    info(0, "\nFailed:", "\n", "\n ".join(failed))

    import pickle

    pickle.dump(
        (eta_bins, psrap, logs),
        open(os.path.splitext(__file__)[0] + ".pkl", "wb"),
        protocol=-1,
    )

    assert passed == set([g.__name__ for g in models])
