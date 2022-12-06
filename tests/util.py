from particle import Particle, ParticleNotFound, InvalidParticle
import typing as _tp
import numpy as np
from multiprocessing import Pool
from os import environ


def reference_charge(pid):
    if isinstance(pid, _tp.Iterable):
        return np.fromiter((reference_charge(pidi) for pidi in pid), np.double)

    try:
        return Particle.from_pdgid(pid).charge
    except (ParticleNotFound, InvalidParticle):
        return np.nan


def run_in_separate_process(fn, *args, timeout=60):
    # Some models need to initialize same fortran code, which can only be
    # initialized once. As a workaround, we run each model in a separate
    # thread. When running several jobs, maxtasksperchild=1 is needed to
    # use a fresh interpreter for each task (not needed here, but still).
    debug = int(environ.get("DEBUG", "0"))
    if debug >= 10:
        out = fn(*args)
    else:
        with Pool(1, maxtasksperchild=1) as p:
            r = p.apply_async(fn, args)
            try:
                out = r.get(timeout=timeout)
            except TimeoutError:
                # usually happens when model aborts and kills child process
                raise TimeoutError("check stdout for errors")

            # In case there any other exception, it can be useful to run in
            # main thread to use the debugger. Set DEBUG=1 to enable this.
            except Exception:
                if debug:
                    out = fn(*args)
                else:
                    raise

    return out
