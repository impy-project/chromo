import os
import numpy as np
from multiprocessing import Pool, freeze_support
import tempfile
import pickle
from impy.constants import GeV
from impy.kinematics import EventKinematics
from impy import impy_config
from impy.util import info
from impy.models import (
    Sibyll23d,
    Sibyll23c,
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


impy_config["user_frame"] = "center-of-mass"


models = [
    Sibyll23d,
    Sibyll23c,
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

xlab_bins = np.linspace(0, 1, 21)


def run_generator(model):

    nevents = 5000

    event_kinematics = EventKinematics(
        ecm=7000 * GeV,
        particle1=2212,
        particle2=2212
        # particle2=(14,7)
    )

    hist_p = 0
    hist_pi = 0
    try:
        log = tempfile.mkstemp()[1]
        generator = model(event_kinematics, logfname=log)
        for event in generator(nevents):
            event.filter_final_state_charged()

            hist_p += np.histogram(
                event.xlab[event.pid == 2212],
                bins=xlab_bins,
                weights=event.xlab[event.pid == 2212] ** 1.7,
            )[0]

            hist_pi += np.histogram(
                event.xlab[np.abs(event.pid) == 211],
                bins=xlab_bins,
                weights=event.xlab[np.abs(event.pid) == 211] ** 1.7,
            )[0]

        return True, model.__class__.__name__, log, hist_p, hist_pi
    except Exception:
        return False, model.__class__.__name__, log, hist_p, hist_pi


def test_all_lab():
    freeze_support()
    pool = Pool(processes=32)
    result = [pool.apply_async(run_generator, (gen,)) for gen in models]
    result = [res.get(timeout=100000) for res in result]

    logs = {}
    xlab_protons = {}
    xlab_piplus = {}
    failed = set()
    passed = set()

    for r, gen, log, hist_p, hist_pi in result:
        if r:
            passed.append(gen)
            xlab_protons[gen] = hist_p
            xlab_piplus[gen] = hist_pi
        else:
            failed.append(gen)

        with open(log) as f:
            logs[gen] = f.read()

    info(0, "Test results for 158 GeV pC collisions in lab frame:\n")
    info(0, "Passed:", "\n", "\n ".join(passed))
    info(0, "\nFailed:", "\n", "\n ".join(failed))

    with open(os.path.splitext(__file__)[0] + ".pkl", "wb") as f:
        pickle.dump((xlab_bins, xlab_protons, xlab_piplus, logs), f, protocol=-1)

    assert passed == set(x.__name__ for x in models)
