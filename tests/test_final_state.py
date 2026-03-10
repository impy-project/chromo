import os
import platform

import boost_histogram as bh
import numpy as np
import pytest
from particle import literals as lp

from chromo import models as im
from chromo.constants import GeV, long_lived
from chromo.kinematics import CenterOfMass
from chromo.util import get_all_models, pdg2name

from .util import run_in_separate_process

# remove omega- from the list of long-lived hadrons
# because the probability is less than 1e3
long_lived_hadrons = [p for p in long_lived if abs(p) not in [13, lp.Omega_minus.pdgid]]


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
    if Model in [im.Sibyll23dStarMixed, im.Sibyll23eStarMixed]:
        pytest.skip(
            reason="SIBYLL* handles decays internally "
            "and ignores most of the decay settings. "
            "Therefore, it doesn't produce short-lived particles."
        )

    if Model is im.Sophia20:
        kin = CenterOfMass(1000 * GeV, "gamma", "p")
    else:
        kin = CenterOfMass(1000 * GeV, "p", "p")
    counts = run_in_separate_process(
        run_model,
        Model,
        kin,
        (
            1000
            if any(
                issubclass(Model, cls) for cls in (im.EposLHC, im.UrQMD34, im.QGSJetIII)
            )
            else 20000
        ),
    )

    # Known issues:
    # SIBYLL-2.1 and UrQMD produce no Omega-
    # QGSJet family produce no Omega-, Xi0, Xi-, (Sigma+, Sigma- only for QGSJetIII)

    known_issues = np.zeros_like(counts, dtype=bool)
    if Model == im.Sibyll21:
        for i, pid in enumerate(long_lived_hadrons):
            if abs(pid) == lp.Omega_minus.pdgid:
                known_issues[i] = True
    elif Model in [im.QGSJet01d, im.QGSJetII03, im.QGSJetII04]:
        for i, pid in enumerate(long_lived_hadrons):
            if abs(pid) in (
                lp.Xi_0.pdgid,
                lp.Xi_minus.pdgid,
                lp.Sigma_plus.pdgid,
                lp.Sigma_minus.pdgid,
            ):
                known_issues[i] = True
    elif Model is im.QGSJetIII:
        for i, pid in enumerate(long_lived_hadrons):
            if abs(pid) in (
                lp.Xi_0.pdgid,
                lp.Xi_minus.pdgid,
            ):
                known_issues[i] = True
    elif issubclass(Model, im.EposLHC):
        for i, pid in enumerate(long_lived_hadrons):
            if pid == lp.Omega_plus_bar.pdgid:
                known_issues[i] = True

    msg = f"{Model.pyname} zero counts for "
    msg += ", ".join(
        f"{pdg2name(pid)}({pid})"
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
