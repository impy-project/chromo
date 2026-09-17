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


def capture_native_printout(model_class, ecm, p1, p2):
    """Generate one event and return the output of print_native_event().

    The generator backends write their printout on the Fortran/C++ level,
    where the output stays buffered unless the process exits, so we run
    the whole thing in a fresh subprocess and capture its stdout. This
    also avoids the restriction that each model can be initialized only
    once per process.
    """
    import pickle
    import subprocess
    import sys

    payload = pickle.dumps((model_class, float(ecm), p1, p2))
    script = """
import pickle, sys
cls, ecm, p1, p2 = pickle.load(sys.stdin.buffer)
from chromo.kinematics import CenterOfMass
generator = cls(CenterOfMass(ecm, p1, p2), seed=1)
for event in generator(1):
    pass
generator.print_native_event()
sys.stdout.flush()
"""
    proc = subprocess.run(
        [sys.executable, "-c", script],
        input=payload,
        capture_output=True,
        timeout=600,
    )
    assert proc.returncode == 0, proc.stderr.decode(errors="replace")[-2000:]
    return proc.stdout.decode(errors="replace")


def baryon_number(pid):
    """Baryon number of a stack record via the particle package.

    Baryons get +-1 via the is_baryon/anti_flag chain, nuclei contribute
    their mass number, everything else zero. Slow; tests only.
    """
    from particle import Particle
    from particle.pdgid import PDGID

    pid = int(pid)
    pg = PDGID(pid)
    if pg == 99999:  # DPMJET hadronization chain placeholder
        return 0
    if pg.is_nucleus:
        return -pg.A if pid < 0 else pg.A
    if pg.is_baryon:
        barred = Particle.from_pdgid(pid).anti_flag.name == "Barred"
        return -1 if (pid < 0 and barred) else 1
    return 0


def charge_number(pdgid):
    """Charge of a stack record in units of e, nuclei handled via pdg2AZ."""
    from particle import Particle

    from chromo.util import pdg2AZ

    pdgid = int(pdgid)
    a, z = pdg2AZ(abs(pdgid))
    if a == 0:
        try:
            return int(Particle.from_pdgid(pdgid).charge)
        except Exception:
            return 0
    return -z if pdgid < 0 else z
