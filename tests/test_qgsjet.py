from functools import lru_cache

import numpy as np
import pytest
from numpy.testing import assert_allclose
from particle import literals as lp

from chromo.common import CrossSectionData
from chromo.constants import GeV
from chromo.kinematics import CenterOfMass
from chromo.models import QGSJet01d, QGSJetII03, QGSJetII04, QGSJetIII

from .util import (
    reference_charge,
    run_in_separate_process,
)

qgsII04_pp_cs_100 = CrossSectionData(
    total=68.21778094810318,
    inelastic=53.97816145792067,
    elastic=14.239619490182513,
    prod=np.nan,
    quasielastic=np.nan,
    coherent=np.nan,
    diffractive_xb=2.331623335930477,
    diffractive_ax=2.0351121873895397,
    diffractive_xx=np.nan,
    diffractive_axb=np.nan,
    diffractive_sum=np.nan,
    b_elastic=16.1990543289839,
)
qgsII03_pp_cs_100 = CrossSectionData(
    total=69.86744771403204,
    inelastic=55.22878210676279,
    elastic=14.638665607269246,
    prod=np.nan,
    quasielastic=np.nan,
    coherent=np.nan,
    diffractive_xb=1.3519578936441186,
    diffractive_ax=1.261095187016783,
    diffractive_xx=np.nan,
    diffractive_axb=np.nan,
    diffractive_sum=np.nan,
    b_elastic=15.989022431287488,
)
qgsIII_pp_cs_100 = CrossSectionData(
    total=65.71870426215214,
    inelastic=52.532563800267894,
    elastic=13.18614046188425,
    prod=np.nan,
    quasielastic=np.nan,
    coherent=np.nan,
    diffractive_xb=2.31753487999351,
    diffractive_ax=1.9929141639968437,
    diffractive_xx=np.nan,
    diffractive_axb=np.nan,
    diffractive_sum=np.nan,
    b_elastic=16.253385724158488,
)
qgs01d_pp_cs_100 = CrossSectionData(
    total=68.03887906922947,
    inelastic=54.24937934684394,
    elastic=13.789499722385532,
    prod=np.nan,
    quasielastic=np.nan,
    coherent=np.nan,
    diffractive_xb=3.099119341661887,
    diffractive_ax=3.099119341661887,
    diffractive_xx=0.6965111778689864,
    diffractive_axb=np.nan,
    diffractive_sum=np.nan,
    b_elastic=np.nan,
)


def run_pp_collision():
    evt_kin = CenterOfMass(100 * GeV, "proton", "proton")
    m = QGSJetII04(evt_kin, seed=4)
    m.set_stable(lp.pi_0.pdgid, True)
    for event in m(1):
        pass
    return event


def run_ab_collision():
    evt_kin = CenterOfMass(100 * GeV, (4, 2), (12, 6))
    m = QGSJetII04(evt_kin, seed=1)
    for event in m(1):
        pass
    return event


@pytest.fixture
@lru_cache(maxsize=1)
def event():
    return run_in_separate_process(run_pp_collision)


@pytest.fixture
@lru_cache(maxsize=1)
def event_ion():
    return run_in_separate_process(run_ab_collision)


def test_impact_parameter(event_ion):
    assert event_ion.impact_parameter > 0


def test_n_wounded(event):
    assert event.n_wounded == (1, 1)


def test_n_wounded_ion(event_ion):
    assert event_ion.n_wounded[0] >= 1
    assert event_ion.n_wounded[1] >= 1


def run_cross_section(Model, p1, p2):
    evt_kin = CenterOfMass(1e3 * GeV, p1, p2)
    m = Model(evt_kin, seed=1)
    return m.cross_section()


@pytest.mark.parametrize(
    "Model, reference_cross_section",
    [
        (QGSJet01d, qgs01d_pp_cs_100),
        (QGSJetII03, qgsII03_pp_cs_100),
        (QGSJetII04, qgsII04_pp_cs_100),
        (QGSJetIII, qgsIII_pp_cs_100),
    ],
)
def test_cross_section(Model, reference_cross_section):
    c = run_in_separate_process(run_cross_section, Model, "p", "p")
    assert c.__eq__(reference_cross_section, rtol=1e-3)


def test_charge(event):
    expected = reference_charge(event.pid)
    ma = np.isnan(expected)
    assert np.mean(ma) < 0.8
    event.charge[ma] = np.nan
    assert_allclose(event.charge, expected)


# def test_vertex(event):
#     assert np.sum(event.vt != 0) > 0


# def test_daughters(event):
#     assert event.daughters.shape == (len(event), 2)
#     # some particles have no daughters
#     assert sum(x[0] == -1 and x[1] == -1 for x in event.daughters) > 0

#     # somes particles have single daughters (elastic scattering in parton shower)
#     assert sum(x[0] >= 0 and x[1] == -1 for x in event.daughters) > 0

#     # some particles have multiple daughters
#     assert sum(x[0] >= 0 and x[1] >= 0 for x in event.daughters) > 0


# def test_mothers(event):
#     assert event.mothers.shape == (len(event), 2)
#     # same particles have no mothers
#     assert sum(x[0] == -1 and x[1] == -1 for x in event.mothers) > 0

#     # most particles have a single mother
#     assert sum(x[0] >= 0 and x[1] == -1 for x in event.mothers) > 0
