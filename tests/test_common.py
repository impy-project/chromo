from chromo.common import CrossSectionData
import numpy as np
import dataclasses


def test_cross_section_empty():
    CrossSectionData()


def test_cross_section_some_filled():
    cx = CrossSectionData(total=12)
    for key, val in dataclasses.asdict(cx).items():
        if key == "total":
            assert val == 12
        else:
            assert np.isnan(val)


def test_cross_section_eq():
    cx1 = CrossSectionData(total=12)
    cx2 = CrossSectionData(total=12)
    assert cx1 == cx2
    cx2.inelastic = 1
    assert cx1 != cx2
    cx2.inelastic = np.nan
    assert cx1 == cx2
    cx2.inelastic = 10
    assert cx1 != cx2
