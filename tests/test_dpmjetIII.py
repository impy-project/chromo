import numpy as np
import chromo
from .util import run_in_separate_process


def run_cross_section_precision():
    event_kin = chromo.kinematics.FixedTarget(1e4, "proton", "O16")
    event_generator = chromo.models.DpmjetIII193(event_kin)

    air = chromo.util.CompositeTarget([("N", 0.78), ("O", 0.22)])
    event_kin = chromo.kinematics.FixedTarget(1e3, "proton", air)

    default_precision = 1000
    # Check the default precision
    assert event_generator.get_cross_section_precision() == default_precision

    # Set a new one
    other_precision = 58
    event_generator.set_cross_section_precision(other_precision)
    assert event_generator.get_cross_section_precision() == other_precision

    # Set to a default precision
    event_generator.set_cross_section_precision()
    assert event_generator.get_cross_section_precision() == default_precision

    trials = 10

    # With small precision
    event_generator.set_cross_section_precision(1)
    cross_section_run1 = np.empty(trials, dtype=np.float64)
    for i in range(trials):
        cross_section_run1[i] = event_generator.cross_section(event_kin).inelastic

    # With default precision
    event_generator.set_cross_section_precision()
    cross_section_run2 = np.empty(trials, dtype=np.float64)
    for i in range(trials):
        cross_section_run2[i] = event_generator.cross_section(event_kin).inelastic

    # Standard deviation should be large for small precision
    assert np.std(cross_section_run1) > np.std(cross_section_run2)


def test_cross_section_precision():
    run_in_separate_process(run_cross_section_precision)
