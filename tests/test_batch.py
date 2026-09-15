import sys

import numpy as np
import pytest

from chromo.kinematics import (
    CenterOfMass,
    CompositeTarget,
    EventFrame,
    EventKinematicsWithRestframe,
    FixedTarget,
    Momentum,
)
from chromo.util import get_all_models

from .util import run_in_separate_process

pytestmark = pytest.mark.skipif(
    sys.platform == "win32", reason="batch models not built on windows"
)


def get_model(pyname):
    return next((c for c in get_all_models() if c.pyname == pyname), None)


def particle_arrays_equal(a, b):
    for field in ("pid", "status", "px", "py", "pz", "en", "m"):
        assert np.array_equal(getattr(a, field), getattr(b, field)), field
    return True


def check_association(events, states):
    # one event per state, in input order, reporting the requested state
    assert len(events) == len(states)
    for event, kin in zip(events, states):
        assert len(event) > 0
        if isinstance(kin.p2, CompositeTarget):
            assert any(event.kin.p2 == c for c in kin.p2.components)
        else:
            assert event.kin == kin


def run_batch_pythia8():
    import chromo.models as im

    states = [
        CenterOfMass(900, "proton", "proton"),
        CenterOfMass(14000, "proton", "proton"),
        CenterOfMass(100, "pi+", "proton"),
        # particle-list style input via beam momenta, cf. issue #188
        EventKinematicsWithRestframe(
            "proton", "proton", beam=(500, -500), frame=EventFrame.CENTER_OF_MASS
        ),
        CenterOfMass(2000, "K+", "neutron"),
    ]
    gen = im.Pythia8(states[1], seed=1)

    events = gen.generate_batch(states, seed=42)
    check_association(events, states)

    # Pythia8 re-seeds its RANMAR from the instance seed on every init,
    # so batch and one-by-one from the same instance must agree exactly
    one_by_one = []
    for kin in states:
        gen.kinematics = kin
        one_by_one.append(next(gen(1)).copy())
    assert len(one_by_one) == len(events)
    for a, b in zip(events, one_by_one):
        particle_arrays_equal(a, b)

    # scrambled stack with duplicates still associates correctly
    dup_stack = [states[i] for i in (4, 2, 1, 0, 2, 0, 3)]
    events4 = gen.generate_batch(dup_stack, seed=7)
    check_association(events4, dup_stack)
    return True


def run_batch_dpmjet():
    import chromo.models as im

    # initialize at the highest energy of the batch, see issue #242
    states = [
        CenterOfMass(50, "proton", "proton"),
        CenterOfMass(300, "proton", "proton"),
        CenterOfMass(1000, "proton", "proton"),
        CenterOfMass(20, "pi+", "proton"),
        CenterOfMass(50, "proton", "neutron"),
        CenterOfMass(300, "proton", "proton"),  # duplicate of states[1]
    ]
    gen = im.DpmjetIII193(states[2], seed=1)

    events = gen.generate_batch(states, seed=42)
    check_association(events, states)

    # duplicates share a single state switch but continue the RNG
    # stream, so their events differ
    assert not np.array_equal(events[1].en, events[5].en)

    # for hadron-proton states, switching the kinematics consumes no
    # random numbers (tabulated cross sections), so a batch whose input
    # is already sorted into groups reproduces one-by-one generation
    # from the same seed exactly
    uniq = sorted(states[:5], key=hash)
    gen.random_state = np.random.default_rng(42).bit_generator.state
    one_by_one = []
    for kin in uniq:
        gen.kinematics = kin
        one_by_one.append(next(gen(1)).copy())
    events2 = gen.generate_batch(uniq, seed=42)
    for a, b in zip(events2, one_by_one):
        particle_arrays_equal(a, b)
    return True


def run_batch_urqmd():
    import chromo.models as im

    states = [
        FixedTarget(Momentum(20), "proton", "proton"),
        FixedTarget(Momentum(100), "proton", "proton"),
        FixedTarget(Momentum(5), "pi+", "proton"),
        FixedTarget(Momentum(20), "proton", "N14"),
    ]
    gen = im.UrQMD34(states[1], seed=3)

    events = gen.generate_batch(states, seed=3)
    assert len(events) == len(states)
    for event, kin in zip(events, states):
        assert len(event) > 0
        assert event.kin.plab == pytest.approx(kin.plab)
        assert event.kin.p1 == kin.p1
        assert event.kin.p2 == kin.p2

    # events are not degenerate
    assert not np.array_equal(events[0].en, events[1].en)
    # note: unlike DPMJET and Pythia8, UrQMD keeps hidden state (nuclear
    # cross section tables, cascade caches) across kinematics switches,
    # so batches are not reproducible across calls with the same seed
    # once a different state was generated in between; only association
    # of events with initial states is guaranteed
    return True


def run_batch_composite_urqmd():
    import chromo.models as im

    air = CompositeTarget((("N", 0.78), ("O", 0.22)), label="air")
    kin0 = FixedTarget(Momentum(100), "proton", "N14")
    states = [FixedTarget(Momentum(100), "proton", air)] * 4
    gen = im.UrQMD34(kin0, seed=5)
    events = gen.generate_batch(states, seed=5)
    assert len(events) == len(states)
    # composite targets are sampled per event, like in generator(nevents)
    allowed = {int(c) for c in air.components}
    for event in events:
        assert int(event.kin.p2) in allowed
    return True


def run_batch_generic_errors(Model):
    good = CenterOfMass(900, "proton", "proton")
    gen = Model(good, seed=1)
    bad_target = CenterOfMass(900, "proton", "Pb208")
    with pytest.raises(ValueError):
        gen.generate_batch([good, bad_target])
    too_low = CenterOfMass(1, "proton", "proton")
    with pytest.raises(ValueError):
        gen.generate_batch([good, too_low])
    # empty batch is fine and returns nothing
    assert gen.generate_batch([]) == []
    return True


def test_batch_pythia8():
    if get_model("Pythia8") is None:
        pytest.skip("Pythia8 not built")
    assert run_in_separate_process(run_batch_pythia8) is True


def test_batch_dpmjet():
    if get_model("DpmjetIII193") is None:
        pytest.skip("DPMJET 19.3 not built")
    assert run_in_separate_process(run_batch_dpmjet) is True


@pytest.mark.parametrize("fn", [run_batch_urqmd, run_batch_composite_urqmd])
def test_batch_urqmd(fn):
    if get_model("UrQMD34") is None:
        pytest.skip("UrQMD 3.4 not built")
    assert run_in_separate_process(fn) is True


def test_batch_generic():
    # batch inherits the normal kinematics validation
    model = get_model("Pythia8")
    if model is None:
        pytest.skip("Pythia8 not built")
    assert run_in_separate_process(run_batch_generic_errors, model) is True