import dataclasses

import numpy as np
import pytest

import chromo.models as im
from chromo.constants import GeV
from chromo.kinematics import (
    CenterOfMass,
    CompositeTarget,
)
from chromo.util import get_all_models

from .util import run_in_separate_process


def run_model(Model, kin):
    try:
        gen = Model(kin, seed=1)
    except ValueError:
        return None

    return (gen.cross_section(), gen._inel_or_prod_cross_section)


@pytest.mark.parametrize("target", ("p", "air"))
@pytest.mark.parametrize("projectile", ("gamma", "pi-", "K+", "p", "He"))
@pytest.mark.parametrize("Model", get_all_models())
def test_generator(projectile, target, Model):
    p1 = projectile
    p2 = target
    if Model is im.UrQMD34:
        pytest.skip("UrQMD34 cross sections broken.")
    if Model is im.Pythia8Angantyr and p1 == "He":
        pytest.skip(
            "Pythia8Angantyr He projectile has no precomputed tables; init too slow for CI"
        )
    if p2 == "air":
        # cannot use Argon in SIBYLL, so make air from N, O only
        p2 = CompositeTarget((("N", 0.78), ("O", 0.22)))

    # increase energy for K+ and SIBYLL23c to avoid prod cs = nan region
    # this is specifig to SIBYLL23C
    if Model == im.Sibyll23c and projectile == "K+" and target == "air":
        kin = CenterOfMass(1e3 * GeV, p1, p2)
    else:
        kin = CenterOfMass(100 * GeV, p1, p2)

    ret = run_in_separate_process(run_model, Model, kin)
    if ret is None:
        assert (
            abs(kin.p1) not in Model.projectiles
            or abs(kin.p2) not in Model.targets
            or kin.ecm < getattr(Model, "_ecm_min", 0)
            or not Model._pair_allowed(abs(kin.p1), abs(kin.p2))
        )
        return
    cs, private_prod = ret
    cstup = dataclasses.astuple(cs)
    assert np.any(cstup), "No valid cross sections returned"
    if p1 == "gamma":
        assert np.nansum(cstup) > 0.1, "Cross section too small" + str(cs)
    else:
        assert np.nansum(cstup) > 20, "Cross section too small" + str(cs)

    assert private_prod is not None and private_prod > 0, "No private cross section set"


# Models that compute the elastic slope and/or rho (Re/Im of the forward
# elastic amplitude) for hadron-nucleon collisions internally
_MODELS_WITH_B_ELASTIC = {
    "DpmjetIII19",
    "EposLHC",
    "EposLHCR",
    "EposLHCRHadrRescattering",
    "Phojet112",
    "Phojet19",
    "Pythia8",
    "QGSJetII03",
    "QGSJetII04",
    "QGSJetIII",
    "Sibyll21",
    "Sibyll23c",
    "Sibyll23d",
    "Sibyll23dStarMixed",
    "Sibyll23e",
    "Sibyll23eStarMixed",
}
_MODELS_WITH_RHO = {
    "Pythia8",
    "Sibyll21",
    "Sibyll23c",
    "Sibyll23d",
    "Sibyll23dStarMixed",
    "Sibyll23e",
    "Sibyll23eStarMixed",
}


@pytest.mark.parametrize("Model", get_all_models())
def test_elastic_amplitude(Model):
    kin = CenterOfMass(1000 * GeV, "p", "p")
    if Model is im.UrQMD34:
        pytest.skip("UrQMD34 cross sections broken.")
    ret = run_in_separate_process(run_model, Model, kin)
    if ret is None:
        assert (
            abs(kin.p1) not in Model.projectiles
            or abs(kin.p2) not in Model.targets
            or not Model._pair_allowed(abs(kin.p1), abs(kin.p2))
        )
        return
    cs, _ = ret
    if Model.pyname in _MODELS_WITH_B_ELASTIC:
        # B in GeV^-2; pp data at 1 TeV: B ~ 16 GeV^-2
        assert 12 < cs.b_elastic < 20, cs
        # optical theorem bound sigma_el <= sigma_tot^2 (1 + rho^2) / (16 pi B)
        # holds with equality for a pure exponential; allow model deviations
        rho = 0 if np.isnan(cs.rho) else cs.rho
        b_mb = cs.b_elastic * 0.389379
        assert cs.elastic == pytest.approx(
            cs.total**2 * (1 + rho**2) / (16 * np.pi * b_mb), rel=0.25
        )
    else:
        assert np.isnan(cs.b_elastic), cs
    if Model.pyname in _MODELS_WITH_RHO:
        # pp data at 1 TeV: rho ~ 0.14
        assert 0.05 < cs.rho < 0.25, cs
    else:
        assert np.isnan(cs.rho), cs


def run_cross_section_list(Model, kin, kins):
    gen = Model(kin, seed=1)
    return [gen.cross_section(k) for k in kins]


@pytest.mark.parametrize(
    "Model",
    [m for m in get_all_models() if m.pyname.startswith("Sibyll")],
)
def test_elastic_amplitude_not_for_nuclei(Model):
    # slope and rho describe the hadron-nucleon amplitude only; they must not
    # be filled (or summed) for nuclear and composite targets
    kin = CenterOfMass(1000 * GeV, "p", "p")
    kins = [
        CenterOfMass(1000 * GeV, "p", "O"),
        CenterOfMass(1000 * GeV, "p", CompositeTarget((("p", 0.5), ("O", 0.5)))),
    ]
    for cs in run_in_separate_process(run_cross_section_list, Model, kin, kins):
        assert np.isnan(cs.b_elastic), cs
        assert np.isnan(cs.rho), cs
        assert cs.prod > 0


def test_non_additive_fields():
    from chromo.common import CrossSectionData

    a = CrossSectionData(total=40.0, b_elastic=15.0, rho=0.1)
    b = CrossSectionData(total=300.0, b_elastic=50.0, rho=0.2)
    s = CrossSectionData()
    s._mul_radd(0.5, a)
    s._mul_radd(0.5, b)
    assert s.total == pytest.approx(170.0)
    assert np.isnan(s.b_elastic)
    assert np.isnan(s.rho)
