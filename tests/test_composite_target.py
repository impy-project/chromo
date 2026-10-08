import numpy as np
import pytest

from chromo.common import CrossSectionData, MCRun
from chromo.kinematics import CenterOfMass, CompositeTarget
from chromo.util import dry_air


def toy_cross_section(p2):
    # p+A toy model: inelastic ~ A^0.7, production slightly below inelastic
    inel = 40.0 * p2.A**0.7
    prod = 0.9 * inel if p2.A > 1 else inel
    return CrossSectionData(inelastic=inel, prod=prod)


class ToyRun(MCRun):
    """Pure-Python generator that only exercises the CompositeTarget logic."""

    _name = "Toy"
    _version = "1"
    _library_name = "toy"
    _event_class = None
    _frame = None
    _ecm_min = 0

    def __init__(self, kin, seed=1):
        self._rng = np.random.default_rng(seed)
        self.n_cross_section_calls = 0
        self.kinematics = kin

    def _generate(self):
        return True

    def _set_kinematics(self, kin):
        pass

    def _cross_section(self, kin=None, max_info=False):
        self.n_cross_section_calls += 1
        return toy_cross_section((kin or self.kinematics).p2)

    def _set_stable(self, pdgid, stable):
        pass

    def plan(self, nevents):
        """Return number of events per component as planned by _composite_plan."""
        result = {}
        for k in self._composite_plan(nevents):
            result[int(self.kinematics.p2)] = k
        return result


COMPONENTS = [("p", 0.5), ("N", 0.3), ("O", 0.2)]


def test_composite_target_weighting_argument():
    t = CompositeTarget(COMPONENTS)
    assert t.weighting == "cross_section"
    t2 = CompositeTarget(COMPONENTS, weighting="number")
    assert t2.weighting == "number"
    assert t != t2
    assert t2.copy() == t2
    assert t2.copy().weighting == "number"
    assert "weighting='number'" in repr(t2)
    assert "weighting" not in repr(t)
    assert dry_air(weighting="number").weighting == "number"
    with pytest.raises(ValueError, match="weighting"):
        CompositeTarget(COMPONENTS, weighting="mass")


def test_composite_target_event_fractions():
    t = CompositeTarget(COMPONENTS)
    sigma = np.array([40.0, 250.0, 280.0])
    expected = t.fractions * sigma / np.sum(t.fractions * sigma)
    assert np.allclose(t.event_fractions(sigma), expected)
    # number fractions are unchanged, they describe the material
    assert np.allclose(t.fractions, [0.5, 0.3, 0.2])

    t2 = CompositeTarget(COMPONENTS, weighting="number")
    assert np.allclose(t2.event_fractions(sigma), t2.fractions)
    assert np.allclose(t2.event_fractions(), t2.fractions)

    with pytest.raises(ValueError):
        t.event_fractions()
    with pytest.raises(ValueError):
        t.event_fractions([1.0, 2.0])
    with pytest.raises(ValueError):
        t.event_fractions([1.0, np.nan, 2.0])
    with pytest.raises(ValueError):
        t.event_fractions([0.0, 0.0, 0.0])


@pytest.mark.parametrize("weighting", ("cross_section", "number"))
def test_composite_cross_section(weighting):
    t = CompositeTarget(COMPONENTS, weighting=weighting)
    gen = ToyRun(CenterOfMass(100, "p", t))
    cs = [toy_cross_section(c) for c in t.components]
    sigma_gen = np.array([cs[0].inelastic, cs[1].prod, cs[2].prod])

    # cross section per target atom, independent of weighting
    xs = gen.cross_section()
    assert xs.inelastic == pytest.approx(sum(t.fractions * [c.inelastic for c in cs]))
    assert gen._inel_or_prod_cross_section == pytest.approx(t.fractions @ sigma_gen)

    if weighting == "number":
        expected = t.fractions
    else:
        expected = t.fractions * sigma_gen / (t.fractions @ sigma_gen)
    assert np.allclose(gen._composite_fractions(gen.kinematics), expected)


@pytest.mark.parametrize("weighting", ("cross_section", "number"))
def test_composite_plan_sampling(weighting):
    t = CompositeTarget(COMPONENTS, weighting=weighting)
    gen = ToyRun(CenterOfMass(100, "p", t), seed=42)
    sigma_gen = np.array(
        [toy_cross_section(t.components[0]).inelastic]
        + [toy_cross_section(c).prod for c in t.components[1:]]
    )
    if weighting == "number":
        expected = t.fractions
    else:
        expected = t.fractions * sigma_gen / np.sum(t.fractions * sigma_gen)

    n = 100_000
    plan = gen.plan(n)
    got = np.array([plan[int(c)] for c in t.components]) / n
    assert np.sum(got) == pytest.approx(1)
    sigma = np.sqrt(expected * (1 - expected) / n)
    assert np.all(np.abs(got - expected) < 5 * sigma)
    # kinematics is restored after the plan
    assert gen.kinematics.p2 == t


def test_composite_fractions_cached():
    gen = ToyRun(CenterOfMass(100, "p", CompositeTarget(COMPONENTS)))
    n = gen.n_cross_section_calls
    gen.plan(10)
    gen.plan(10)
    assert gen.n_cross_section_calls == n


def test_composite_fractions_fallback_on_invalid_cross_section():
    class NanRun(ToyRun):
        def _cross_section(self, kin=None, max_info=False):
            return CrossSectionData()

    t = CompositeTarget(COMPONENTS)
    with pytest.warns(RuntimeWarning, match="number fractions"):
        gen = NanRun(CenterOfMass(100, "p", t))
    assert np.allclose(gen._composite_fractions(gen.kinematics), t.fractions)
    assert np.isnan(gen._inel_or_prod_cross_section)


def test_composite_fractions_follow_kinematics_change():
    gen = ToyRun(CenterOfMass(100, "p", CompositeTarget(COMPONENTS)))
    f_xs = gen._composite_fractions(gen.kinematics)
    assert not np.allclose(f_xs, [0.5, 0.3, 0.2])
    gen.kinematics = CenterOfMass(
        100, "p", CompositeTarget(COMPONENTS, weighting="number")
    )
    assert np.allclose(gen._composite_fractions(gen.kinematics), [0.5, 0.3, 0.2])
