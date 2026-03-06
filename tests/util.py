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
    import os
    import sys

    # On Windows, Fortran's WRITE(*,...) writes to OS-level fd 1. If pytest
    # redirected fd 1 to a pipe (for output capture), the pipe buffer can fill
    # while the parent waits for the subprocess, causing a deadlock. Redirect
    # fd 1 to devnull in the subprocess to prevent this.
    if sys.platform == "win32":
        _devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(_devnull, 1)
        os.close(_devnull)

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
        # Read from the queue FIRST (with timeout), then wait for the process.
        # On Windows (and generally), the subprocess's Queue Feeder thread can
        # block if the parent doesn't drain the queue while waiting for the
        # process to die — a classic pipe-buffer deadlock. Reading first avoids
        # this.
        try:
            out = queue.get(timeout=timeout)
        except Exception:
            p.kill()
            assert False, "subprocess timed out or queue get failed"
        p.join(timeout=5)
        if p.is_alive():
            p.kill()
        p.close()
    return out
