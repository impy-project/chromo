from pathlib import Path

import numpy as np
import pytest
import uproot
import yaml
from numpy.testing import assert_allclose, assert_equal

from chromo.common import CrossSectionData, EventData
from chromo.kinematics import CompositeTarget, EventKinematicsWithRestframe
from chromo.writer import Root


def make_event(n):
    assert n >= 2
    rng = np.random.default_rng(1)
    pid = rng.choice([211, 130, 2212], size=n)
    pid[:2] = 2212
    status = np.arange(n)
    mothers = rng.choice(n, size=(n, 2))
    mothers[:2] = 0
    mothers[:, 1] = 0
    return EventData(
        ("foo", "1.0"),
        EventKinematicsWithRestframe("p", "He", beam=(-3, 4)),
        1,
        1.0,
        (2, 3),
        100.0,
        pid,
        status,
        1.1 + status,
        2.2 + status,
        3.3 + status,
        4.4 + status,
        5.5 + status,
        6.6 + status,
        7.7 + status,
        8.8 + status,
        9.9 + status,
        10.01 + status,
        mothers,
        None,
    )


@pytest.mark.parametrize("write_vertices", (False, True))
@pytest.mark.parametrize("overflow", (False, True))
@pytest.mark.parametrize("target", ("He", "air"))
def test_Root(write_vertices, overflow, target):
    if target == "air":
        air = CompositeTarget([("N", 0.75), ("O", 0.25)])
        kin = EventKinematicsWithRestframe("p", air, plab=5)
    else:
        kin = EventKinematicsWithRestframe("p", target, beam=(-3, 4))

    class Model:
        label: str = "foo"
        seed = 1
        kinematics = kin

        def cross_section(self):
            return CrossSectionData(total=6.6)

    events = [
        make_event(2),
        make_event(5),
        make_event(4),
    ]
    model = Model()

    # name must contain all parameters to not cause collisions when test is run parallel
    p = Path(f"test_writer_{write_vertices}_{overflow}_{target}.root")

    writer = Root(p, model, write_vertices=write_vertices, buffer_size=5)
    if overflow:
        with writer:
            for event in events:
                writer.write(event)
            with pytest.raises(RuntimeError):
                writer.write(make_event(10))
    else:
        with writer:
            for event in events:
                writer.write(event)

    with uproot.open(p) as f:
        tree = f["event"]
        data = yaml.safe_load(tree.title)
        ref = {
            "seed": 1,
            "projectile_id": 2212,
            "projectile_momentum": -3 if target == "He" else 5,
            "target_id": 1000020040 if target == "He" else repr(air),
            "target_momentum": 4 if target == "He" else 0,
            "model": "foo",
            "sigma_total": 6.6,
            "energy_unit": "GeV",
            "sigma_unit": "mb",
        }
        if write_vertices:
            ref["length_unit"] = "mm"
        assert data == ref
        assert tree.num_entries == len(events)
        d = tree.arrays()
        if not write_vertices:
            assert "vx" not in tree
        for i, event in enumerate(events):
            assert_equal(d["pdgid"][i], event.pid[2:])
            assert_allclose(d["px"][i], event.px[2:])
            assert_equal(d["parent"][i], np.maximum(event.mothers[2:, 0] - 2, -1))
            if write_vertices:
                assert_allclose(d["vx"][i], event.vx[2:])

    p.unlink()
