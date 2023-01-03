from impy.kinematics import CenterOfMass
from impy.remote_control import MCRunRemote as MCRun
from impy.common import CrossSectionData
from impy.util import EventFrame
import pytest


class DummyEvent:
    def __init__(self, generator, kinematics, counter=0):
        self.kinematics = kinematics
        if generator is None:
            self.counter = counter
        else:
            self.counter = generator.event.counter

    def copy(self):
        return DummyEvent(None, self.kinematics, self.counter)

    def __repr__(self):
        return f"DummyEvent(None, {self.kinematics}, {self.counter})"

    def __eq__(self, other):
        return self.kinematics == other.kinematics and self.counter == other.counter


class Dummy(MCRun):
    _name = "Dummy"
    _version = "1.0"
    _event_class = DummyEvent
    _library_name = "_eposlhc"
    _frame = EventFrame.CENTER_OF_MASS

    def _once(self, raise_at=None):
        self.event = DummyEvent(None, None)
        self.raise_at = raise_at

    def _generate(self):
        if self.raise_at is not None:
            if self.raise_at == 0:
                raise ValueError("from raise_at")
            self.raise_at -= 1
        self.event.counter += 1
        # fake drawing a random number
        self._lib.crranma4.ntot += 1
        return True

    def _cross_section(self, kin):
        self._lib.crranma4.ntot += 1
        return CrossSectionData(inelastic=float(self._lib.crranma4.ntot))

    def _set_kinematics(self, kin):
        self._kin = kin

    def _set_stable(self, pid, stable):
        pass

    def __repr__(self):
        return f"<Dummy at 0x{id(self):x}>"


@pytest.mark.parametrize("timeout", (100, 0))
def test_dummy_1(timeout):
    model = Dummy(seed=1, timeout=timeout)
    assert model.random_state.counter == 0

    kin = CenterOfMass(100, "p", "p")
    events = []
    for event in model(kin, 5):
        if timeout:
            assert model.random_state.counter == 0
        else:
            assert model.random_state.counter > 0
        events.append(event)

    expected = [DummyEvent(None, kin, i + 1) for i in range(5)]
    assert events == expected

    assert model.random_state.counter == 5

    c = model.cross_section(kin)
    assert c.inelastic == 5


@pytest.mark.parametrize("timeout", (100, 0))
def test_dummy_2(timeout):
    model = Dummy(seed=1, timeout=timeout)

    kin = CenterOfMass(100, "p", "p")
    events = []
    for i, event in enumerate(model(kin, 10)):
        if timeout:
            assert model.random_state.counter == 0
        else:
            assert model.random_state.counter > 0
        if i == 3:
            break
        events.append(event)

    expected = [DummyEvent(None, kin, i + 1) for i in range(3)]
    assert events == expected

    assert model.random_state.counter >= 3


@pytest.mark.parametrize("timeout", (100, 0))
def test_dummy_3(timeout):
    model = Dummy(seed=1, timeout=timeout, raise_at=3)

    kin = CenterOfMass(100, "p", "p")
    events = []

    with pytest.raises(ValueError):
        for event in model(kin, 10):
            if timeout:
                assert model.random_state.counter == 0
            else:
                assert model.random_state.counter > 0
            events.append(event)

    expected = [DummyEvent(None, kin, i + 1) for i in range(3)]
    assert events == expected

    assert model.random_state.counter >= 3


@pytest.mark.parametrize("timeout", (100, 0))
def test_dummy_4(timeout):
    model = Dummy(seed=1, timeout=timeout)

    kin = CenterOfMass(100, "p", "p")
    events = []

    with pytest.raises(ValueError):
        for i, event in enumerate(model(kin, 10)):
            if timeout:
                assert model.random_state.counter == 0
            else:
                assert model.random_state.counter > 0
            if i == 3:
                raise ValueError()
            events.append(event)

    expected = [DummyEvent(None, kin, i + 1) for i in range(3)]
    assert events == expected

    assert model.random_state.counter >= 3
