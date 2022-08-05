from __future__ import print_function

import sys
import os
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(root_dir)

from impy.definitions import *
from impy.constants import *
from impy.kinematics import EventKinematics
from impy import impy_config, pdata
from impy.util import info
from impy.writer import HepMCWriter
import pytest

event_kinematics = EventKinematics(ecm=7000 * GeV, p1pdg=2212, p2pdg=2212)

impy_config["user_frame"] = "laboratory"


@pytest.mark.parametrize(
    "model_tag",
    [
        "SIBYLL23D",
        # "DPMJETIII306",
        # "EPOSLHC"
    ],
)
def test_hepmc_writer(model_tag):
    # To run this test do `pytest tests/test_hepmc_writer.py`
    # This test fails because the event record written by HepMC3 C++ is bad,
    # a lot of particles are missing. Either a bug in the original impy record or a bug in the
    # HepMC3 C++ code (not the pyhepmc code).
    generator = make_generator_instance(interaction_model_by_tag[model_tag])
    generator.init_generator(event_kinematics)

    test_file = "test_hepmc_writer_file_%s.dat" % model_tag

    event_data = []
    with HepMCWriter(test_file) as w:
        for event in generator.event_generator(event_kinematics, 3):
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

    import pyhepmc

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
    import os

    os.unlink(test_file)
