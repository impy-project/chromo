from chromo.common import CrossSectionData, EventData, MCEvent
from chromo.kinematics import CenterOfMass, EventFrame
import numpy as np
import dataclasses
import pickle
from types import SimpleNamespace
import pytest
from numpy.testing import assert_equal


@pytest.fixture
def evt():
    i = np.array([1, 2, 3])
    f = np.array([1.1, 2.2, 3.3])
    p = np.array([[1, -1], [2, -1], [2, 3]])
    return EventData(
        "generator",
        CenterOfMass(10, "p", "p"),
        1,
        0.5,
        (1, 1),
        i,
        i,
        f,
        f,
        f,
        f,
        f,
        f,
        f,
        f,
        f,
        f,
        p,
        p,
    )


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


def test_EventData_copy_and_pickle(evt):
    evt2 = evt.copy()
    assert evt2 == evt

    s = pickle.dumps(evt)
    evt3 = pickle.loads(s)

    assert evt3 == evt


class DummyEvent(MCEvent):
    def __init__(self):
        hepevt = SimpleNamespace(
            nevhep=1,
            nhep=2,
            idhep=np.ones(3, dtype=np.int32),
            isthep=np.ones(3, dtype=np.int32),
            phep=np.ones((3, 5), dtype=np.double).T,
            vhep=np.ones((3, 4), np.double).T,
            jmohep=np.ones((3, 2), np.int32),
            jdahep=np.ones((3, 2), np.int32),
        )

        lib = SimpleNamespace(hepevt=hepevt)

        generator = SimpleNamespace(
            _lib=lib,
            name="foo",
            version="bar",
            kinematics=CenterOfMass(10, "p", "p"),
            _frame=EventFrame.CENTER_OF_MASS,
        )

        super().__init__(generator)

    def _charge_init(self, npart):
        return np.zeros(npart)

    def _add_init_beam_info(self):
        self._prepend_initial_beam()


def test_MCEvent_copy():
    evt = DummyEvent()
    evt2 = evt.copy()
    assert evt2 == evt

    s = pickle.dumps(evt)
    evt3 = pickle.loads(s)

    assert evt3 == evt


def test_EventData_select(evt):
    x = evt[1]
    assert x.pid == evt.pid[1]

    x = evt[[True, False, True]]
    assert_equal(x.pid, [1, 3])
