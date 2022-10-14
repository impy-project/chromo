from particle import Particle, ParticleNotFound, InvalidParticle
import typing as _tp
import numpy as np
from multiprocessing import Pool
import os
from impy import models as im
from impy.common import MCRun
import pytest


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
    with Pool(1, maxtasksperchild=1) as p:
        r = p.apply_async(fn, args)
        try:
            out = r.get(timeout=timeout)
        except TimeoutError:
            # usually happens when model aborts and kills child process
            raise TimeoutError("check stdout for errors")

        # In case there any other exception, it can be useful to run in
        # main thread to get proper backtrace. Uncomment the following to
        # enable this.

        # except Exception:
        #     fn(*args)

    return out


# remove this when git lfs issue is fixed
def skip_on_ci_if_model_is_incompatible(Model):
    if os.environ.get("CI", False) and Model in (
        im.QGSJet01d,
        im.QGSJetII03,
        im.QGSJetII04,
        im.Phojet191,
        im.Phojet112,
        im.Phojet193,
        im.EposLHC,
        im.DpmjetIII306,
        im.DpmjetIII191,
        im.DpmjetIII193,
    ):
        pytest.skip(
            "model cannot succeed on CI, because git lfs does not work",
            allow_module_level=True,
        )


def get_all_models(module):
    result = []
    for key in dir(module):
        obj = getattr(module, key)
        try:
            # fails if obj is not a class
            if issubclass(obj, MCRun):
                result.append(obj)
        except TypeError:
            pass
    return result
