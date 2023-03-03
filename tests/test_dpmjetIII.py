import numpy as np
import chromo
from .util import run_in_separate_process


def run_cross_section_ntrials():
    event_kin = chromo.kinematics.FixedTarget(1e4, "proton", "O16")
    event_generator = chromo.models.DpmjetIII193(event_kin)

    air = chromo.util.CompositeTarget([("N", 0.78), ("O", 0.22)])
    event_kin = chromo.kinematics.FixedTarget(1e3, "proton", air)

    default_precision = 1000
    # Check the default precision
    assert event_generator.glauber_trials == default_precision

    # Set a new one
    other_precision = 58
    event_generator.glauber_trials = other_precision
    assert event_generator.glauber_trials == other_precision

    trials = 10

    # With small precision
    event_generator.glauber_trials = 1
    cross_section_run1 = np.empty(trials, dtype=np.float64)
    for i in range(trials):
        cross_section_run1[i] = event_generator.cross_section(event_kin).inelastic

    # With default precision
    event_generator.glauber_trials = 1000
    cross_section_run2 = np.empty(trials, dtype=np.float64)
    for i in range(trials):
        cross_section_run2[i] = event_generator.cross_section(event_kin).inelastic

    # Standard deviation should be large for small precision
    assert np.std(cross_section_run1) > np.std(cross_section_run2)


def test_cross_section_ntrials():
    run_in_separate_process(run_cross_section_ntrials)
