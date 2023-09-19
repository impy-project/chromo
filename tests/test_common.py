from chromo.common import CrossSectionData, EventData, MCEvent
from chromo.kinematics import CenterOfMass, EventFrame
import numpy as np
import dataclasses
import pickle
from types import SimpleNamespace
import pytest
from contextlib import nullcontext
from numpy.testing import assert_equal
from chromo.util import get_all_models


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


def test_non_diffractive():
    csd = CrossSectionData(inelastic=10.0, diffractive_xb=3.0, diffractive_sum=8.0)
    assert csd.non_diffractive == 2.0


def test_diffractive():
    csd = CrossSectionData(
        diffractive_xb=3.0, diffractive_ax=2.0, diffractive_xx=1.0, diffractive_axb=4.0
    )
    assert csd.diffractive == 10.0

    csd.diffractive_sum = 5.0
    assert csd.diffractive == 5.0


def test_mul_radd():
    csd1 = CrossSectionData(total=5.0, inelastic=10.0)
    csd2 = CrossSectionData(total=10.0, inelastic=5.0)

    csd1._mul_radd(2, csd2)

    assert_equal(csd1.total, 25.0)
    assert_equal(csd1.inelastic, 20.0)


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


@pytest.mark.parametrize("Model", get_all_models())
def test_models_beam(Model):
    """Tests whether all models have correct beam particles"""
    evt_kin = CenterOfMass(100, "proton", "proton")
    if Model.pyname in "Sophia20":
        evt_kin = CenterOfMass(100, "photon", "proton")
    elif Model.name in ["DPMJET-III", "EPOS"]:
        evt_kin = CenterOfMass(100, "N", "O")
    elif Model.name in ["SIBYLL"]:
        evt_kin = CenterOfMass(100, "p", "O")

    with pytest.warns(RuntimeWarning) if Model.name in [
        "DPMJET-III",
        "UrQMD",
    ] else nullcontext():
        generator = Model(evt_kin, seed=1)
        for event in generator(100):
            event.kin._beam_data = None
            beam = event.kin._get_beam_data(event.kin.frame)
            event.kin._beam_data = None
            for field, beam_field in beam.items():
                event_field = getattr(event, field)
                if event_field is None or field == "daughters":
                    continue
                assert np.allclose(
                    event_field[0:2], beam_field
                ), f"{field}: {np.allclose(event_field[0:2], beam_field)}, {event_field[0:2]}, {beam_field}"
