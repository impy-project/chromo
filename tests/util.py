from particle import Particle, ParticleNotFound, InvalidParticle
import typing as _tp
import numpy as np
from os import environ
import time


def reference_charge(pid):
    if isinstance(pid, _tp.Iterable):
        return np.fromiter((reference_charge(pidi) for pidi in pid), np.double)

    try:
        return Particle.from_pdgid(pid).charge
    except (ParticleNotFound, InvalidParticle):
        return np.nan


def _target(queue, fn, args):
    out = fn(*args)
    queue.put(out)


def run_in_separate_process(fn, *args, timeout=120):
    import multiprocessing as mp

    # Some models need to initialize same fortran code, which can only be
    # initialized once. As a workaround, we run each model in a separate
    # thread. When running several jobs, maxtasksperchild=1 is needed to
    # use a fresh interpreter for each task (not needed here, but still).
    debug = int(environ.get("DEBUG", "0"))
    if debug >= 10:
        out = fn(*args)
    else:
        queue = mp.Queue()
        p = mp.Process(target=_target, args=(queue, fn, args))
        p.start()
        for _ in range(timeout * 2):
            if p.is_alive():
                time.sleep(0.5)
            else:
                break
        if queue.empty():
            assert False, "queue empty, process probably crashed"
        out = queue.get(timeout=1)
    return out
