import numpy as np
from multiprocessing import Pool, freeze_support, TimeoutError
import tempfile
import pickle
from pathlib import PurePath as Path
from impy.constants import GeV
from impy.kinematics import EventKinematics
from impy import impy_config
import impy.models
import abc


def run_generator(model, nevents):
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
    norm = 1.0 / float(nevents) / eta_widths
    hist = np.zeros(len(eta_centers))
    try:
        log = tempfile.mkstemp()[1]
        generator = model(event_kinematics)
        for event in generator(nevents):
            event.filter_final_state_charged()
            hist += np.histogram(event.eta, bins=eta_bins)[0]
        return True, model.__name__, log, eta_bins, hist * norm
    except Exception:
        return False, model.__name__, log, eta_bins, hist * norm


def test_all_cms():
    freeze_support()

    # put all models in a list
    models = []
    for key in dir(impy.models):
        obj = getattr(impy.models, key)
        if isinstance(obj, abc.ABCMeta):
            models.append(obj)

    result = []
    with Pool(processes=len(models)) as pool:
        jobs = [pool.apply_async(run_generator, (model, 100)) for model in models]
        for job in jobs:
            try:
                r = job.get(timeout=1000)
                result.append(r)
            except TimeoutError:
                pass

    failed = set([x.__name__ for x in models])
    passed = set()
    psrap = {}
    logs = {}

    eta_bins = result[0][3]

    for r, gen, log, _, hist in result:
        if r:
            passed.add(gen)
            psrap[gen] = hist

        with open(log) as f:
            logs[gen] = f.read()

    failed -= passed

    with open(Path(__file__).stem + ".pkl", "wb") as f:
        pickle.dump((eta_bins, psrap, logs), f)

    assert not failed
