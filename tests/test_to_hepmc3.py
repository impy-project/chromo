from collections import defaultdict
from impy.kinematics import CenterOfMass
from impy import models as im
from impy.constants import TeV
from .util import run_in_separate_process, xfail_on_ci_if_model_is_incompatible
import numpy as np
import pytest


def make_event(Model):
    ekin = CenterOfMass(1 * TeV, 2212, 2212)
    m = Model(ekin, seed=1)
    for event in m(100):
        if len(event) > 10:  # to skip elastic events
            break
    return event  # MCEvent is pickeable, but restored as EventData


@pytest.mark.parametrize("Model", (im.Sibyll21, im.Pythia6, im.EposLHC))
def test_to_hepmc3(Model):
    xfail_on_ci_if_model_is_incompatible(Model)

    event = run_in_separate_process(make_event, Model)

    unique_vertices = defaultdict(list)
    for i, pa in enumerate(event.parents):
        assert pa.shape == (2,)
        if np.all(pa == 0):
            continue
        pa = (pa[0], pa[1])
        # normalize intervals
        if pa[1] == 0:
            pa = (pa[0], pa[0])
        unique_vertices[pa].append(i)

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

    unique_vertices2 = defaultdict(list)
    for v in hev.vertices:
        pi = [p.id for p in v.particles_in]
        if len(pi) == 1:
            pa = (pi[0], pi[0])
        else:
            pa = (min(pi), max(pi))
        children = [p.id - 1 for p in v.particles_out]
        unique_vertices2[pa] = children

    assert unique_vertices == unique_vertices2
