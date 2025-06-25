import time
import typing as _tp
from os import environ

import numpy as np
from particle import InvalidParticle, Particle, ParticleNotFound


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


def run_in_separate_process(fn, *args, timeout=600):
    import multiprocessing as mp

    # Some models need to initialize same fortran code, which can only be
    # initialized once. As a workaround, we run each model in a separate
    # Process.
    debug = abs(int(environ.get("DEBUG", "0")))
    if debug >= 10:
        out = fn(*args)
    else:
        ctx = mp.get_context("spawn")
        queue = ctx.Queue()
        p = ctx.Process(target=_target, args=(queue, fn, args))
        p.start()
        step = 0.5
        for _ in range(int(timeout / step)):
            if p.is_alive():
                time.sleep(step)
            else:
                break
        if queue.empty():
            assert False, "queue empty, process probably crashed"
        out = queue.get(timeout=1)
        p.join()
        p.close()
    return out
