<<<<<<< HEAD:tests/test_all_cms.py
import numpy as np
from multiprocessing import Pool, freeze_support, TimeoutError
import tempfile
import pickle
from pathlib import PurePath as Path
from impy.constants import GeV
from impy.kinematics import EventKinematics
from impy import impy_config
from impy.common import MCRun
import impy.models
import traceback


def run_generator(model, nevents):

    if model.__name__ == "Sophia20":
        event_kinematics = EventKinematics(elab=7000 * GeV, p1pdg=22, p2pdg=2212)
        impy_config["user_frame"] = "center-of-mass"
    else:
        
    if (model.__name__== "Sophia20"):        
        event_kinematics = EventKinematics(
            elab=7000 * GeV,
            p1pdg=22,
            p2pdg=2212
        )
        impy_config["user_frame"] = "center-of-mass"
    else:
        event_kinematics = EventKinematics(
                ecm=7000 * GeV,
                p1pdg=2212,
                p2pdg=2212
                # nuc2_prop=(14,7)
            )

            impy_config["user_frame"] = "center-of-mass"

    model_name = model.__name__

    eta_bins = np.linspace(-5, 5, 21)
    eta_widths = eta_bins[1:] - eta_bins[:-1]
    eta_centers = 0.5 * (eta_bins[1:] + eta_bins[:-1])
    norm = 1.0 / float(nevents) / eta_widths
    hist = np.zeros(len(eta_centers))
    success = False
    with tempfile.NamedTemporaryFile(mode="w+") as log:
        try:
            log.write(model_name + " init\n")
            generator = model(event_kinematics, logfname=log.name)
            log.write(model_name + " start\n")
            for event in generator(nevents):
                event.filter_final_state_charged()
                hist += np.histogram(event.eta, bins=eta_bins)[0]
            log.write(model_name + " finish\n")
            success = True
        except Exception:
            log.write(model_name + " abort\n")
            log.write(traceback.format_exc())

        log.flush()
        log.seek(0)
        return success, model_name, log.read(), eta_bins, hist * norm


def test_all_cms():
    freeze_support()

    # put all models in a list
    models = []
    for key in dir(impy.models):
        obj = getattr(impy.models, key)
        if hasattr(obj, "__mro__") and MCRun in obj.__mro__:
            models.append(obj)

    result = []
    with Pool(1, maxtasksperchild = 1) as pool:
        jobs = [pool.apply_async(run_generator, (model, 100)) for model in models]
        for job in jobs:
            try:
                r = job.get(timeout=30)
                result.append(r)
            except TimeoutError:
                pass

    failed = set([x.__name__ for x in models])
    psrap = {}
    logs = {}

    for r, model_name, log, eta_bins, hist in result:
        if r:
            failed.remove(model_name)
            psrap[model_name] = hist
        else:
            logs[model_name] = log
            print(log)

    with open(Path(__file__).stem + ".pkl", "wb") as f:
        pickle.dump((eta_bins, psrap, logs), f)

    assert failed == set()
>>>>>>> 914028d (All tests in test_all_cms.py pass):impy/tests/test_all_cms.py
