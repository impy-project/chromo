"""
This example script generates phi mesons for proton-proton collisions at 13 TeV.

You can specify the model and how many phi mesons you want generated on the command line.
Beware that some models (the QGSJet series) provide no particle history, which means that
this model produces no phis at all. Running the script with a QGSJet model will stall
forever.

The script generates a ROOT file with three trees.
    phi:
        Contains the eta and pt / MeV of phis, whether the phi was prompt (it actually
        always is), and the number of charged particles in the event in which the phi
        was produced.
    run_info:
        Contains only a single number per branch, the information on how many inelastic
        events had to be generated to produce that many phis and the inelastic
        cross-section. With this information, one can compute the differential
        cross-section for phi production for the model.
    charged:
        Contains the distributed of charged particles for all generated events.
"""

from chromo.kinematics import CenterOfMass
from chromo.constants import MeV, TeV, long_lived
from chromo.util import get_all_models
from tqdm import tqdm
import numpy as np
import numba as nb
from particle import literals as lp
import uproot
import argparse


parser = argparse.ArgumentParser()
parser.add_argument(
    "model_spec", help="specify model to run using substring of its pyname"
)
parser.add_argument("number", type=int, help="number of phis to generate")

args = parser.parse_args()


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
    if args.model_spec.lower() in Model.pyname.lower():
        matches.append(Model)
assert len(matches) == 1


def run(Model, requested):
    kin = CenterOfMass(13 * TeV, "p", "p")
    phi_pid = lp.phi_1020.pdgid

    m = Model(kin, seed=1)

    eta = []
    pt = []
    prompt = []
    charged = []
    all_charged = []
    n_events = 0
    with tqdm(total=requested, desc=m.pyname, ascii=True) as bar:
        while True:
            for event in m(100):
                n_events += 1
                n_phi = np.sum(event.pid == phi_pid)
                nch = len(event.final_state_charged())
                all_charged.append(nch)
                if n_phi > 0:
                    if np.all(event.pid[:2] == 2212) and np.all(
                        event.mothers[:2, 0] == -1
                    ):
                        event = event[2:]  # cut beam particles
                    ma = event.pid == phi_pid
                    imot = event.mothers[:, 0]
                    eta += list(event.eta[ma])
                    pt += list(event.pt[ma] / MeV)
                    prompt += [
                        is_prompt(idx, event.pid, imot)
                        for idx in np.arange(len(event))[ma]
                    ]
                    charged += [nch for _ in range(n_phi)]
                    bar.update(n_phi)
                    if len(eta) >= requested:
                        with uproot.recreate(f"phi_{m.pyname}.root") as f:
                            f["phi"] = {
                                "eta": eta,
                                "pt": pt,
                                "prompt": prompt,
                                "charged": charged,
                            }
                            f["run_info"] = {
                                "n_events": np.array([n_events]),
                                "inelastic_cross_section": np.array(
                                    [m.cross_section().inelastic]
                                ),
                            }
                            f["charged"] = {"charged": all_charged}
                            return


run(matches[0], args.number)
