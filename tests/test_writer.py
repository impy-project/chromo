from impy.writer import Root
from impy.common import CrossSectionData, EventData
from impy.kinematics import CenterOfMass
from types import SimpleNamespace
import numpy as np
import uproot
from numpy.testing import assert_equal
import yaml


def test_Root():
    config = SimpleNamespace(
        seed=1,
        projectile_id=2,
        projectile_momentum=3.3,
        target_id=4,
        target_momentum=5.5,
        model=SimpleNamespace(label="foo"),
    )
    c = CrossSectionData(total=6.6)

    ia = np.array([211, 130, 2212])
    fa = np.array([1.5, 2.4, 3.5])
    i2a = np.array([[1, 2], [0, 0], [2, 1]])

    event = EventData(
        ("foo", "1.0"),
        CenterOfMass(10, "p", "p"),
        "user-frame",
        1,
        1.0,
        (2, 3),
        ia,
        ia,
        fa,
        fa,
        fa,
        fa,
        fa,
        fa,
        fa,
        fa,
        fa,
        fa,
        i2a,
        None,
    )

    with Root("foo.root", config, c) as f:
        f.write(event)

    with uproot.open("foo.root") as f:
        tree = f["event"]
        data = yaml.safe_load(tree.title)
        assert data == {
            "seed": 1,
            "projectile_id": 2,
            "projectile_momentum": 3.3,
            "target_id": 4,
            "target_momentum": 5.5,
            "model": "foo",
            "sigma_total": 6.6,
            "energy_unit": "MeV",
            "length_unit": "mm",
            "sigma_unit": "mb",
        }
        assert tree.num_entries == 1
        assert_equal(tree["pdgid"], ia[2:])
