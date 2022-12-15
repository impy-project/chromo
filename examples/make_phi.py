from impy import models as im
from impy.kinematics import CenterOfMass
from impy.constants import MeV, TeV, long_lived
from impy.util import get_all_models
from tqdm import tqdm
import joblib
import numpy as np
import numba as nb
from particle import literals as lp
import uproot
import sys

model_spec = sys.argv[1]
requested = int(sys.argv[2])


@nb.vectorize
def is_long_lived(pid):
    """
    Return true if PDGID is a long-lived particle.

    Parameters
    ----------
    pid: int or array-like
        PDGID of particle.

    Returns
    -------
    bool or array of bool

    Notes
    -----
    This function can be used in numba-jitted functions and in normal Python code.
    When used from Python, it is vectorized. Inside a numba-jitted function, it
    accepts a single pid.
    """
    return abs(pid) in long_lived


@nb.njit
def is_prompt(idx, pid, imot):
    """
    Return true for prompt particles.

    Parameters
    ----------
    idx : int
        Index of current particle.
    pid : array
        Array containing the PDGID of particles.
    imot : array
        Array containing the index of the parent particles or -1
        if there is no parent.

    Returns
    -------
    bool
    """
    if len(pid) != len(imot):
        raise ValueError("pid and imot do not have same length")

    while idx != -1:
        idx = imot[idx]
        if idx == -1:
            return True
        if idx >= len(pid):
            raise IndexError("idx out of bounds")
        if is_long_lived(pid[idx]):
            return False

    return True


matches = []
for Model in get_all_models():
    if model_spec.lower() in Model.pyname.lower():
        matches.append(Model)
assert len(matches) == 1


def run(Model, requested):
    kin = CenterOfMass(13 * TeV, "p", "p")
    phi_pid = lp.phi_1020.pdgid

    m = Model(kin, seed=1)

    eta = []
    pt = []
    prompt = []
    n_events = 0
    with tqdm(total=requested, desc=m.pyname, ascii=True) as bar:
        while True:
            for event in m(100):
                n_events += 1
                n_phi = np.sum(event.pid == phi_pid)
                if n_phi > 0:
                    event = event[2:]  # cut beam particles
                    ma = event.pid == phi_pid
                    imot = event.parents[:, 0] - 1
                    eta += list(event.eta[ma])
                    pt += list(event.pt[ma] / MeV)
                    prompt += [
                        is_prompt(idx, event.pid, imot)
                        for idx in np.arange(len(event))[ma]
                    ]
                    bar.update(n_phi)
                    if len(eta) >= requested:
                        with uproot.recreate(f"phi_{m.pyname}.root") as f:
                            f["phi"] = {"eta": eta, "pt": pt, "prompt": prompt}
                            f["run_info"] = {
                                "n_events": np.array([n_events]),
                                "inelastic_cross_section": np.array(
                                    [m.cross_section().inelastic]
                                ),
                            }
                            return


run(matches[0], requested)
