import numpy as np
import chromo
from .util import run_in_separate_process


def run_cross_section_ntrials():
    event_kin = chromo.kinematics.FixedTarget(1e4, "proton", "O16")
    event_generator = chromo.models.DpmjetIII193(event_kin)

    air = chromo.util.CompositeTarget([("N", 0.78), ("O", 0.22)])
    event_kin = chromo.kinematics.FixedTarget(1e3, "proton", air)

    default_ntrials = 1000
    # Check the default ntrials
    assert event_generator.get_hA_AA_glauber_trials() == default_ntrials

    # Set a new one
    other_ntrials = 58
    event_generator.set_hA_AA_glauber_trials(other_ntrials)
    assert event_generator.get_hA_AA_glauber_trials() == other_ntrials

    # Set to a default ntrials
    event_generator.set_hA_AA_glauber_trials()
    assert event_generator.get_hA_AA_glauber_trials() == default_ntrials

    nruns = 10

    # With small ntrials
    event_generator.set_hA_AA_glauber_trials(1)
    cross_section_run1 = np.empty(nruns, dtype=np.float64)
    for i in range(nruns):
        cross_section_run1[i] = event_generator.cross_section(event_kin).inelastic

    # With default ntrials
    event_generator.set_hA_AA_glauber_trials()
    cross_section_run2 = np.empty(nruns, dtype=np.float64)
    for i in range(nruns):
        cross_section_run2[i] = event_generator.cross_section(event_kin).inelastic

    # Standard deviation should be large for small ntrials
    assert np.std(cross_section_run1) > np.std(cross_section_run2)


def test_cross_section_ntrials():
    run_in_separate_process(run_cross_section_ntrials)
