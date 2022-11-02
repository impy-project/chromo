from impy.kinematics import CenterOfMass, FixedTarget
from impy import models as im
from impy.constants import GeV
from impy.models.sophia import Sophia20
from .util import (
    run_in_separate_process,
    get_all_models,
)
import numpy as np
import pytest

# generate list of all models in impy.models
models = get_all_models(im)


def run(Model):
    evt_kin = CenterOfMass(10 * GeV, "proton", "proton")
    if Model is Sophia20:
        evt_kin = CenterOfMass(10 * GeV, "photon", "proton")
    m = Model(evt_kin, seed=1)
    for event in m(100):
        if len(event) > 10:  # to skip small events
            break
    return event  # MCEvent is pickeable, but restored as EventData


@pytest.mark.parametrize("Model", models)
def test_to_hepmc3(Model):

    if Model == im.UrQMD34:
        pytest.xfail("UrQMD34 FAILS, should be FIXED!!!")

    event = run_in_separate_process(run, Model)

    # special case for models that only have final-state particles
    if Model is im.UrQMD34 or Model.name in ("PhoJet", "DPMJET-III"):
        hev = event.to_hepmc3()
        # only final state is stored
        fs = event.final_state()
        assert len(hev.particles) == len(fs)
        for i, p in enumerate(hev.particles):
            assert p.pid == fs.pid[i]
            assert p.status == fs.status[i]
        assert len(hev.vertices) == 0
        return  # test ends here
    # special case for Pythia8, which does not contain the parton show
    elif Model is im.Pythia8:
        # parton shower is skipped
        from impy.constants import quarks_and_diquarks_and_gluons

        ma = True
        apid = np.abs(event.pid)
        for p in quarks_and_diquarks_and_gluons:
            ma &= apid != p
        event = event[ma]

    unique_vertices = {}
    for i, pa in enumerate(event.parents):
        assert pa.shape == (2,)
        if np.all(pa == 0):
            continue
        # normalize intervals
        if pa[1] == 0:
            pa = (pa[0], pa[0])
        pa = (pa[0] - 1, pa[1])
        # in case of overlapping ranges of incoming particles
        # the earlier vertex keeps them
        for (a, b) in unique_vertices:
            if pa != (a, b) and a <= pa[0] < b:
                pa = b, pa[1]
        unique_vertices.setdefault(pa, []).append(i)

    # check that parent ranges do not exceed particle range;
    # that's a requirement for a valid particle history
    nmax = len(event.px)
    for i, (a, b) in enumerate(unique_vertices):
        assert a >= 0 or a == -1
        assert (
            b <= nmax
        ), f"vertex {i} has parent range {(a, b)} which exceeds particle record nmax={nmax}"

    # not all vertices have locations different from zero,
    # create unique fake vertex locations for testing
    for ch in unique_vertices.values():
        i = ch[0]
        event.vx[i] = i
        event.vy[i] = i + 1
        event.vz[i] = i + 2
        event.vt[i] = i + 3

    hev = event.to_hepmc3()

    assert len(hev.run_info.tools) == 1
    assert hev.run_info.tools[0] == (*event.generator, "")
    assert len(hev.particles) == len(event)
    assert len(hev.vertices) == len(unique_vertices)

    for i, p in enumerate(hev.particles):
        assert p.momentum.x == event.px[i]
        assert p.momentum.y == event.py[i]
        assert p.momentum.z == event.pz[i]
        assert p.momentum.e == event.en[i]
        assert p.status == event.status[i]
        assert p.pid == event.pid[i]
        assert p.id == i + 1

    for i, v in enumerate(hev.vertices):
        k = v.particles_out[0].id - 1
        assert v.position.x == event.vx[k]
        assert v.position.y == event.vy[k]
        assert v.position.z == event.vz[k]
        assert v.position.t == event.vt[k]

    hev_vertices = {}
    for v in hev.vertices:
        pi = [p.id - 1 for p in v.particles_in]
        if len(pi) == 1:
            pa = (pi[0], pi[0] + 1)
        else:
            pa = (min(pi), max(pi) + 1)
        children = [p.id - 1 for p in v.particles_out]
        hev_vertices[pa] = children

    assert unique_vertices == hev_vertices
