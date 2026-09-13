import numpy as np
import pytest
from particle import literals as lp

from chromo.constants import GeV
from chromo.kinematics import CenterOfMass
from chromo.models import UrQMD34

from .util import run_in_separate_process


def count_multiplicities(p1, n_events=200):
    """Return final-state multiplicities of n_events collisions of p1 on proton.

    UrQMD can only be initialized once per process, so this is meant to be
    run via run_in_separate_process.
    """
    kin = CenterOfMass(100 * GeV, p1, int(lp.proton.pdgid))
    gen = UrQMD34(kin, seed=1234)
    mults = []
    for event in gen(n_events):
        mults.append(len(event.final_state()))
    return mults


# UrQMD does not generate events with elastic scattering only. Such events
# produced many final states with only the two beam particles for pion
# projectiles, see https://github.com/impy-project/chromo/issues/45
@pytest.mark.parametrize(
    "projectile,pdgid",
    [("pi+", int(lp.pi_plus.pdgid)), ("p", int(lp.proton.pdgid))],
)
def test_no_empty_events_with_baryonic_projectile(projectile, pdgid):
    mults = run_in_separate_process(count_multiplicities, pdgid, 200, timeout=1200)
    mults = np.array(mults)
    assert len(mults) == 200
    low_mult_fraction = np.mean(mults <= 2)
    assert low_mult_fraction < 0.05, (
        f"{projectile}+p at 100 GeV: {100 * low_mult_fraction:.1f}% of events "
        "have <=2 final state particles"
    )


def test_pion_low_mult_comparable_to_proton():
    """Low-multiplicity fraction for pi+p must be comparable to p+p."""
    pi_mults = run_in_separate_process(
        count_multiplicities, int(lp.pi_plus.pdgid), 200, timeout=1200
    )
    p_mults = run_in_separate_process(
        count_multiplicities, int(lp.proton.pdgid), 200, timeout=1200
    )
    f_pi = np.mean(np.array(pi_mults) <= 2)
    f_p = np.mean(np.array(p_mults) <= 2)
    # both are inelastic-only samples now; allow a loose binomial margin
    assert f_pi - f_p < 0.05, f"pi+p low-mult fraction {f_pi} vs p+p {f_p}"
