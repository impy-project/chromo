import os
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from impy.constants import GeV
from impy.kinematics import EventKinematics
from impy.writer import HepMCWriter
from impy.models import Sibyll23d, DpmjetIII306, EposLHC  # noqa
import pytest
import pyhepmc

# TODO for Hans: this test fails because new pyhepmc lacks
# pyhepmc.fill_genevent_from_hepevent


@pytest.mark.parametrize(
    "model",
    [
        Sibyll23d,
        DpmjetIII306,
        EposLHC,
    ],
)
def test_hepmc_writer(model):
    # To run this test do `pytest tests/test_hepmc_writer.py`
    # This test fails because the event record written by HepMC3 C++ is bad,
    # a lot of particles are missing. Either a bug in the original impy record or a
    # bug in the HepMC3 C++ code (not the pyhepmc code).

    ekin = EventKinematics(ecm=7000 * GeV, p1pdg=2212, p2pdg=2212)
    gen = model(ekin)

    test_file = f"test_hepmc_writer_file_{model.__name__}.dat"

    event_data = []
    with HepMCWriter(test_file) as w:
        for event in gen(3):
            # event.filter_final_state()
            n = event.npart
            pem = event._pem_arr.T
            vt = event._vt_arr.T
            pid = event.p_ids
            status = event.status
            parents = {}
            if event.parents is not None:
                for i, (a, b) in enumerate(event.parents.T):
                    par = set()
                    for j in range(a, b):
                        par.add((pid[j], pem[j, 3]))
                    if par:
                        parents[(pid[i], pem[i, 3])] = vt[i], par
            event_data.append(
                (pem.copy(), vt.copy(), pid.copy(), status.copy(), parents.copy())
            )
            w.write(event)

    for ievent, event in enumerate(pyhepmc.open(test_file)):
        assert event is not None
        assert event.event_number == ievent

        pem_ref, vt_ref, pid_ref, status_ref, parents_ref = event_data[ievent]

        particles = event.particles
        n = len(particles)
        pem = np.empty((n, 5))
        pid = np.empty(n, dtype=int)
        status = np.empty(n, dtype=int)
        parents = {}
        for i, p in enumerate(particles):
            assert i + 1 == p.id
            pem[i] = (
                p.momentum.px,
                p.momentum.py,
                p.momentum.pz,
                p.momentum.e,
                p.generated_mass,
            )
            pid[i] = p.pid
            status[i] = p.status
            par = set()
            for q in p.parents:
                par.add((q.pid, q.momentum.e))
            if par:
                parents[(p.pid, p.momentum.e)] = p.production_vertex.position, par

        assert_allclose(pem_ref, pem)
        assert_array_equal(pid_ref, pid)
        assert_array_equal(status_ref, status)
        assert parents == parents_ref

    # delete test_file if test is successful
    os.unlink(test_file)
