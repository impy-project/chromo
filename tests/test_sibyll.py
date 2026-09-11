import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

from chromo.common import CrossSectionData
from chromo.constants import TeV
from chromo.kinematics import CenterOfMass
from chromo.util import get_all_models

from .util import reference_charge, run_in_separate_process

cs_sibyll21 = CrossSectionData(
    total=117.90274047851562,
    inelastic=85.13753509521484,
    elastic=32.76520538330078,
    prod=np.nan,
    quasielastic=np.nan,
    coherent=np.nan,
    diffractive_xb=6.198373,
    diffractive_ax=6.198373,
    diffractive_xx=1.8042221,
    diffractive_axb=0,
    diffractive_sum=np.nan,
    b_elastic=np.nan,
)

cs_sibyll23 = CrossSectionData(
    total=105.7269966641395,
    inelastic=76.78422435244065,
    elastic=28.94277231169886,
    prod=np.nan,
    quasielastic=np.nan,
    coherent=np.nan,
    diffractive_xb=5.985736090273257,
    diffractive_ax=5.985736090273257,
    diffractive_xx=1.662186514630767,
    diffractive_axb=0,
    diffractive_sum=np.nan,
    b_elastic=np.nan,
)


def get_sibylls():
    """Get the list of all Sibylls"""
    return [cl for cl in get_all_models() if cl.pyname.startswith("Sibyll")]


def event_run(model):
    evt_kin = CenterOfMass(10 * TeV, "p", "p")
    m = model(evt_kin, seed=1)  # Use the passed model
    for event in m(1):
        pass
    return event.copy()  # copy is pickleable


@pytest.mark.parametrize("model", get_sibylls())
def test_charge(model):
    event = run_in_separate_process(event_run, model)
    expected = reference_charge(event.pid)
    assert_allclose(event.charge, expected)


@pytest.mark.parametrize("model", get_sibylls())
def test_daughters(model):
    event = run_in_separate_process(event_run, model)
    assert event.daughters is None


# @pytest.mark.parametrize("model", get_sibylls())
# def test_mothers(model):
#     # This stuff doesn't work for SIBYLL and should have never worked
#     event = run_in_separate_process(event_run, model)
#     # check that there are particles with a single parent
#     # and that parent is short-lived or an initial particle

#     ma = (event.mothers[:, 0] > 0) & (event.mothers[:, 1] == -1)

#     assert np.sum(ma) > 0

#     idx = event.mothers[ma, 0]

#     # remove beam particles (if any)
#     ma = event.status[idx] != 3
#     idx = idx[ma]

#     # remaining particles should be short-lived or
#     # bookkeeping particles like K0
#     for pid in event.pid[idx]:
#         p = Particle.from_pdgid(pid)
#         assert p.ctau is None or p.ctau < 1


@pytest.mark.parametrize("model", get_sibylls())
def test_vertex(model):
    event = run_in_separate_process(event_run, model)
    # no vertex info in Sibyll21
    assert_equal(event.vx, 0)
    assert_equal(event.vy, 0)
    assert_equal(event.vz, 0)
    assert_equal(event.vt, 0)


def run_cross_section(model, p1, p2):
    evt_kin = CenterOfMass(10 * TeV, p1, p2)
    m = model(evt_kin, seed=1)
    return m.cross_section()


@pytest.mark.parametrize("model", get_sibylls())
def test_cross_section(model):
    c = run_in_separate_process(run_cross_section, model, "p", "p")
    if model.pyname == "Sibyll21":
        reference_cs = cs_sibyll21
    else:
        reference_cs = cs_sibyll23
    assert c.__eq__(reference_cs, rtol=1e-3)
    assert c.diffractive_axb == 0
    assert_allclose(
        c.non_diffractive,
        c.inelastic - c.diffractive_xb - c.diffractive_ax - c.diffractive_xx,
    )


def run_with_runtime_warning(model, p1, p2):
    evt_kin = CenterOfMass(10 * TeV, p1, p2)
    with pytest.warns(RuntimeWarning) as record:
        _ = model(evt_kin, seed=1)
    return [str(w.message) for w in record]


@pytest.mark.parametrize(
    "model",
    [
        (
            pytest.param(
                m,
                marks=pytest.mark.skip(
                    reason="Sibyll defines all cs for its projectiles"
                ),
            )
            if m.__name__ == "Sibyll21"
            else m
        )
        for m in get_sibylls()
    ],
)
def test_not_supported_cross_section(model):
    msgs = run_in_separate_process(run_with_runtime_warning, model, "D_0", "p")

    assert len(msgs) == 1
    assert "Cross section for <PDGID: 421> projectiles not supported" in msgs[0]


def wounded_check(model, kins, n=20):
    """Run one generator sequentially through a list of kinematics,
    collecting (n_wounded[0], n_wounded[1], impact_parameter) per event.

    The first events of each kinematics are dropped: SIBYLL only fills
    the geometry blocks on non-coherent-diffraction events, so early
    values may be uninitialized. Must run in a separate process
    (Fortran global state).
    """
    out = []
    gen = None
    for kin in kins:
        if gen is None:
            gen = model(kin)
        else:
            gen.kinematics = kin
        gen.set_stable(111)
        events = list(gen(n + 2))
        out.append(
            [
                (int(e.n_wounded[0]), int(e.n_wounded[1]), e.impact_parameter)
                for e in events[2:]
            ]
        )
    return out


@pytest.mark.parametrize("model", get_sibylls())
def test_wounded_hadron_proton(model):
    """h + N has no Glauber geometry: no impact parameter, but one
    wounded nucleon on each side (same convention as other generators)."""
    (results,) = run_in_separate_process(
        wounded_check, model, [CenterOfMass(5 * TeV, "p", "p")]
    )
    for na, nb, b in results:
        assert (na, nb) == (1, 1)
        assert b == 0.0  # SIBYLL does not sample b for hadron-nucleon


@pytest.mark.parametrize("model", get_sibylls())
def test_wounded_proton_nucleus(model):
    """p + A runs through the sibyll() path with the Glauber geometry in
    /S_CNCM0/: projectile side has exactly one wounded nucleon, target
    side fluctuates, impact parameter is sampled."""
    (results,) = run_in_separate_process(
        wounded_check, model, [CenterOfMass(5 * TeV, "p", (16, 8))]
    )
    assert all(na == 1 for na, nb, b in results)
    assert all(1 <= nb <= 16 for na, nb, b in results)
    assert any(nb > 1 for na, nb, b in results), "no multiple interactions"
    assert all(b > 0 for na, nb, b in results)


@pytest.mark.parametrize("model", get_sibylls())
def test_wounded_nucleus_proton(model):
    """O + p runs through the sibnuc() path with geometry in /CNUCMS/.
    chromo reports wounded as (projectile, target): the proton target is
    singly wounded or, for elastic events, not wounded at all."""
    (results,) = run_in_separate_process(
        wounded_check, model, [CenterOfMass(5 * TeV, (16, 8), "p")]
    )
    assert all(0 <= na <= 16 for na, nb, b in results)
    assert all(nb in (0, 1) for na, nb, b in results)
    assert any(na > 1 for na, nb, b in results), "no multiple interactions"
    assert all(b > 0 for na, nb, b in results)


@pytest.mark.parametrize("model", get_sibylls())
def test_wounded_nucleus_nucleus(model):
    """O + O runs through the sibnuc() path with geometry in /CNUCMS/."""
    (results,) = run_in_separate_process(
        wounded_check, model, [CenterOfMass(5 * TeV, (16, 8), (16, 8))]
    )
    assert all(0 <= na <= 16 and 0 <= nb <= 16 for na, nb, b in results)
    assert any(na > 1 and nb > 1 for na, nb, b in results)
    assert all(b > 0 for na, nb, b in results)


@pytest.mark.parametrize("model", get_sibylls())
def test_wounded_no_stale_blocks(model):
    """SIBYLL never clears its geometry common blocks. Retargeting a
    generator must never leak wounded counts or impact parameters from a
    previous collision type into events of the current one."""
    results = run_in_separate_process(
        wounded_check,
        model,
        [
            CenterOfMass(5 * TeV, "p", (16, 8)),  # fill /S_CNCM0/
            CenterOfMass(5 * TeV, "p", "p"),  # must not leak into h + N
            CenterOfMass(5 * TeV, (16, 8), (16, 8)),  # fill /CNUCMS/
            CenterOfMass(5 * TeV, "p", (16, 8)),  # hA geometry is live
        ],
    )
    # after a p + O run, hadron-nucleon events must report exactly
    # (1, 1) and b = 0, not the Glauber values from the previous stage
    for na, nb, b in results[1]:
        assert (na, nb, b) == (1, 1, 0.0)
    # after switching to O + O, the /CNUCMS/ geometry must be freshly sampled
    assert all(b > 0 for na, nb, b in results[2])
    # and switching back to p + O must sample /S_CNCM0/ again
    assert all(na == 1 and nb >= 1 and b > 0 for na, nb, b in results[3])
