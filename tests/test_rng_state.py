import importlib
import pickle
from copy import deepcopy

import numpy as np
import pytest

import chromo.models as im
from chromo.common import MCRun
from chromo.constants import GeV, TeV
from chromo.kinematics import CenterOfMass, FixedTarget
from chromo.util import get_all_models

from .util import run_in_separate_process


def run_rng_state(Model):
    if Model is im.Sophia20:
        evt_kin = FixedTarget(13 * TeV, "photon", "proton")
    elif Model is im.UrQMD34:
        evt_kin = CenterOfMass(50 * GeV, "proton", "proton")
    elif Model in (im.Pythia8Cascade, im.Pythia8Angantyr):
        evt_kin = CenterOfMass(13 * TeV, "proton", "N14")
    elif Model is im.Fluka:
        # Fluka p+p event generation is unstable above ~200 GeV CMS, so
        # use 100 GeV here instead of 13 TeV.
        evt_kin = CenterOfMass(100 * GeV, "proton", "proton")
    else:
        evt_kin = CenterOfMass(13 * TeV, "proton", "proton")

    generator = Model(evt_kin, seed=1)
    if Model is im.Fluka:
        # Pythia8 decay handler's RNG is independent from FLUKA's Ranmar
        # state and would diverge across the two passes; disable it for
        # this reproducibility check.
        generator._activate_decay_handler(on=False)

    nevents = 10
    # EposLHC needs more events to fail
    if Model is im.EposLHC:
        nevents = 50

    # Save a initial state to a variable:
    state_0 = deepcopy(generator.random_state)
    # Generate nevents events
    states = [deepcopy(generator.random_state) for _ in generator(nevents)]

    # Pickle generator state after nevents
    pickled_state_1 = pickle.dumps(generator.random_state)

    # Restore initial state from variable
    generator.random_state = state_0

    # Compare states after each generated event
    for i, _ in enumerate(generator(nevents)):
        assert (
            states[i] == generator.random_state
        ), f"states differ after {i+1} generation"

    state_1a = generator.random_state
    state_1 = pickle.loads(pickled_state_1)
    # check that pickled state_1 and reproduced state_1a are equal
    assert state_1a == state_1


def _event_physics_arrays(event):
    return [
        getattr(event, name)
        for name in (
            "pid",
            "status",
            "charge",
            "px",
            "py",
            "pz",
            "en",
            "m",
            "vx",
            "vy",
            "vz",
            "vt",
        )
    ]


def run_rng_state_events_reproducible(Model):
    """After restoring the initial RNG state, the very same events must be
    produced, not just matching RNG states. For UrQMD34 this also checks
    that the one-time initialization of its embedded Pythia6 is fully
    absorbed by the deterministic warm-up in the constructor (issue #64).
    """
    assert Model is im.UrQMD34
    evt_kin = CenterOfMass(50 * GeV, "proton", "proton")
    generator = Model(evt_kin, seed=1)

    nevents = 10
    state_0 = deepcopy(generator.random_state)
    # deepcopy because events are views into the Fortran event record
    # for some models, and the record is overwritten by the next event
    events = [deepcopy(event) for event in generator(nevents)]

    generator.random_state = state_0
    for i, event in enumerate(generator(nevents)):
        assert np.isclose(
            event.impact_parameter, events[i].impact_parameter
        ), f"impact parameter differs at {i}"
        for a, b in zip(_event_physics_arrays(event), _event_physics_arrays(events[i])):
            assert np.array_equal(a, b), f"events differ at {i}"


def test_rng_state_events_reproducible():
    run_in_separate_process(run_rng_state_events_reproducible, im.UrQMD34)


@pytest.mark.skipif(not hasattr(im, "Pythia6"), reason="Pythia6 not built")
def test_rng_state_pythia6_mpi_averages():
    # With PCG64(101), event 10 updates the running MPI averages of PYMULT
    # enough to change event 4 on replay unless they are part of the state.
    run_in_separate_process(run_rng_state_with_bitgen, im.Pythia6, np.random.PCG64, 101)


def run_rng_state_with_bitgen(Model, bitgen_class, seed):
    """Test RNG state save/restore with specific bit generator."""
    if Model is im.Sophia20:
        evt_kin = FixedTarget(13 * TeV, "photon", "proton")
    elif Model is im.UrQMD34:
        evt_kin = CenterOfMass(50 * GeV, "proton", "proton")
    elif Model in (im.Pythia8Cascade, im.Pythia8Angantyr):
        evt_kin = CenterOfMass(13 * TeV, "proton", "N14")
    elif Model is im.Fluka:
        # See run_rng_state — same Fluka caveats apply.
        evt_kin = CenterOfMass(100 * GeV, "proton", "proton")
    else:
        evt_kin = CenterOfMass(13 * TeV, "proton", "proton")

    rng = np.random.Generator(bitgen_class(seed=seed))
    generator = Model(evt_kin, seed=rng)
    if Model is im.Fluka:
        # See run_rng_state.
        generator._activate_decay_handler(on=False)

    nevents = 10 if Model is not im.EposLHC else 50

    state_0 = deepcopy(generator.random_state)
    states = [deepcopy(generator.random_state) for _ in generator(nevents)]

    generator.random_state = state_0

    for i, _ in enumerate(generator(nevents)):
        current_state = generator.random_state
        # Use str() for comparison to handle numpy arrays in nested dicts
        assert str(states[i]) == str(
            current_state
        ), f"states differ after {i+1} generation with {bitgen_class.__name__}"


@pytest.mark.parametrize("Model", get_all_models())
def test_rng_state(Model):
    if Model in (im.EposLHCR, im.EposLHCRHadrRescattering):
        pytest.skip(
            f"{Model.pyname} maintains UrQMD internal state that affects event "
            "generation but is not captured by RNG state."
        )

    run_in_separate_process(run_rng_state, Model)


@pytest.mark.parametrize("Model", get_all_models())
@pytest.mark.parametrize(
    "bitgen_class,seed",
    [
        (np.random.PCG64, 11111),
        (np.random.MT19937, 22222),
        (np.random.Philox, 33333),
        (np.random.SFC64, 44444),
        (np.random.PCG64DXSM, 55555),
    ],
)
def test_rng_state_bitgens(Model, bitgen_class, seed):
    """Test different NumPy bit generators with all models."""
    if Model in (im.EposLHCR, im.EposLHCRHadrRescattering):
        pytest.skip(
            f"{Model.pyname} maintains UrQMD internal state that affects event "
            "generation but is not captured by RNG state."
        )

    run_in_separate_process(run_rng_state_with_bitgen, Model, bitgen_class, seed)


def run_bitgen_save_restore_roundtrip(bitgen_class):
    # EPOS saves and restores the RNG state through ranfgt/ranfst,
    # which copy the NumPy bit generator state in C (rangen.c).
    lib = importlib.import_module("chromo.models._eposlhc")

    def draw(n):
        out = np.zeros(max(n, 1))
        if n:
            lib.rmmard(out, n, 0)
        return out[:n]

    # vary the position inside the buffered output (Philox) and key (MT19937)
    for nskip in range(9):
        rng = np.random.Generator(bitgen_class(seed=1234))
        lib.npy.bitgen = rng.bit_generator.ctypes.bit_generator.value
        lib.npy.gen_id = MCRun._get_bitgen_id(rng.bit_generator)
        draw(nskip)
        lib.ranfgt(0.0)
        expected = draw(6)
        lib.ranfst(0.0)
        assert np.array_equal(draw(6), expected), f"differs after {nskip} draws"


@pytest.mark.skipif(not hasattr(im, "EposLHC"), reason="EPOS-LHC not built")
@pytest.mark.parametrize(
    "bitgen_class",
    [
        np.random.PCG64,
        np.random.PCG64DXSM,
        np.random.MT19937,
        np.random.Philox,
        np.random.SFC64,
    ],
)
def test_bitgen_save_restore_roundtrip(bitgen_class):
    run_in_separate_process(run_bitgen_save_restore_roundtrip, bitgen_class)
