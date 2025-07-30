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
    pytest.xfail(reason="The SIBYLL charge test needs a fix.")
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


def check_wounded(model, kin):
    gen = model(kin)
    gen.set_stable(111)


    wounded_na_list = []
    wounded_nb_list = []
    wounded_b_list = []
    cnucms_na_list = []
    cnucms_nb_list = []
    cnucms_b_list = []
    s_cncm0_na_list = []
    s_cncm0_b_list = []

    for event in gen(10):

        # Check wounded candidates [Main check]
        wounded_na, wounded_nb = event.n_wounded[0], event.n_wounded[1]
        wounded_b = event.impact_parameter
        wounded_na_list.append(wounded_na !=0)
        wounded_nb_list.append(wounded_nb !=0)
        wounded_b_list.append(wounded_b !=0)

        # Check cnucms status [Detail check]
        cnucms_na, cnucms_nb = event._lib.cnucms.na, event._lib.cnucms.nb
        cnucms_b = event.impact_parameter
        cnucms_na_list.append(cnucms_na !=0)
        cnucms_nb_list.append(cnucms_nb !=0)
        cnucms_b_list.append(cnucms_b !=0)

        # Check s_cncm0 status [Detail check]
        s_cncm0_na_list.append(event._lib.s_cncm0.na !=0)
        s_cncm0_b_list.append(event._lib.s_cncm0.b !=0)

    main_errors = []
    detail_errors = []

    if not any(wounded_na_list):
        main_errors.append(f"{gen.label=} | wounded_na_list is empty")
    if not any(wounded_nb_list):
        main_errors.append(f"{gen.label=} | wounded_nb_list is empty")
    if not any(wounded_b_list):
        main_errors.append(f"{gen.label=} | wounded_b_list is empty")

    if not any(cnucms_na_list):
        detail_errors.append(f"{gen.label=} | cnucms_na_list is empty")
    if not any(cnucms_nb_list):
        detail_errors.append(f"{gen.label=} | cnucms_nb_list is empty")
    if not any(cnucms_b_list):
        detail_errors.append(f"{gen.label=} | cnucms_b_list is empty")
    if not any(s_cncm0_na_list):
        detail_errors.append(f"{gen.label=} | s_cncm0_na_list is empty")
    if not any(s_cncm0_b_list):
        detail_errors.append(f"{gen.label=} | s_cncm0_b_list is empty")

    return main_errors, detail_errors


@pytest.mark.parametrize("model", get_sibylls())
def test_wounded_proton_proton(model):
    pytest.xfail(reason="No impact parameter in proton-proton collision, neither in " \
    "event._lib.s_cncm0.b or event._lib.cnucms.b")
    main_errors, detail_errors = run_in_separate_process(check_wounded, model, CenterOfMass(5 * TeV, "p", "p"))
    assert not main_errors

@pytest.mark.parametrize("model", get_sibylls())
def test_wounded_proton_nucleus(model):
    main_errors, detail_errors = run_in_separate_process(check_wounded, model, CenterOfMass(5 * TeV, "p", (16, 8)))
    assert not main_errors

@pytest.mark.parametrize("model", get_sibylls())
def test_wounded_nucleus_proton(model):
    main_errors, detail_errors = run_in_separate_process(check_wounded, model, CenterOfMass(5 * TeV, (16, 8), "p"))
    assert not main_errors

@pytest.mark.parametrize("model", get_sibylls())
def test_wounded_nucleus_nucleus(model):
    main_errors, detail_errors = run_in_separate_process(check_wounded, model, CenterOfMass(5 * TeV, (16, 8), (16, 8)))
    assert not main_errors
