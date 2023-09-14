from chromo.constants import GeV, long_lived
from chromo.kinematics import CenterOfMass
from chromo import models as im
import pytest
from .util import run_in_separate_process
from chromo.util import get_all_models, pdg2name
import boost_histogram as bh
import numpy as np
import platform
import os
from particle import literals as lp

long_lived_hadrons = [p for p in long_lived if abs(p) != 13]


def run_model(Model, kin, number=1000):
    gen = Model(kin, seed=1)
    # Decay of particles produced by QGSJet01d
    # produce Sigma-, Xi-, Xi+, Xi0
    if Model.pyname.startswith("QGSJet"):
        gen._activate_decay_handler(on=False)
    h = bh.Histogram(bh.axis.IntCategory(long_lived_hadrons))
    for event in gen(number):
        ev = event.final_state()
        h.fill(ev.pid)
    return h.values().astype(int)


@pytest.mark.skipif(
    "CI" in os.environ and platform.system() == "Windows",
    reason="skip to speed up CI on Windows",
)
@pytest.mark.parametrize("Model", get_all_models())
def test_generator(Model):
    if Model == im.Sibyll23StarMixed:
        pytest.skip(
            reason="SIBYLL* handles decays internally "
            "and ignores most of the decay settings. "
            "Therefore, it doesn't produce short-lived particles."
        )

    if Model == im.EposLHC:
        pytest.xfail(
            reason="SHOULD BE FIXED: EposLHC don't to produce"
            " some of the required particles "
            " (the type of particles and frequency of fails"
            " seem depend on debug code)"
            " This should not happen because we use the same seed."
            ""
        )

    if Model is im.Sophia20:
        kin = CenterOfMass(1000 * GeV, "gamma", "p")
    else:
        kin = CenterOfMass(1000 * GeV, "p", "p")
    counts = run_in_separate_process(
        run_model, Model, kin, 1000 if Model in (im.EposLHC, im.UrQMD34) else 20000
    )

    # Known issues:
    # SIBYLL-2.1 and UrQMD produce no Omega-
    # QGSJet family produce no Omega-, Xi0, Xi-, Sigma+, Sigma-
    known_issues = np.zeros_like(counts, dtype=bool)
    if Model == im.Sibyll21:
        for i, pid in enumerate(long_lived_hadrons):
            if abs(pid) == lp.Omega_minus.pdgid:
                known_issues[i] = True
    elif Model == im.UrQMD34:
        for i, pid in enumerate(long_lived_hadrons):
            if pid == lp.Omega_minus.pdgid:
                known_issues[i] = True
    elif Model.pyname.startswith("QGSJet"):
        for i, pid in enumerate(long_lived_hadrons):
            if abs(pid) in (
                lp.Omega_minus.pdgid,
                lp.Xi_0.pdgid,
                lp.Xi_minus.pdgid,
                lp.Sigma_plus.pdgid,
                lp.Sigma_minus.pdgid,
            ):
                known_issues[i] = True

    msg = f"{Model.pyname} zero counts for "
    msg += ", ".join(
        pdg2name(pid)
        for pid, n, ki in zip(long_lived_hadrons, counts, known_issues)
        if n == 0 and not ki
    )
    test = np.all((counts > 0) ^ known_issues)
    assert test, (
        msg
        + ", ".join(
            f"{pid} {n} {ki}"
            for pid, n, ki in zip(long_lived_hadrons, counts, known_issues)
        )
        + f"__: {(counts > 0) ^ known_issues}"
    )
