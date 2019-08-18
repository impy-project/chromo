from __future__ import print_function

import sys
import os
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))
sys.path.append(root_dir)

from impy.definitions import *
from impy.constants import *
from impy.kinematics import EventKinematics
from impy.common import impy_config, pdata
from impy.util import info
from impy.writer import HepMCWriter
import pytest

event_kinematics = EventKinematics(ecm=7000 * GeV,
                                   p1pdg=2212,
                                   p2pdg=2212)

impy_config["user_frame"] = 'laboratory'


@pytest.mark.parametrize("model_tag", ["SIBYLL23C", "DPMJETIII306"])
def test_hepmc_writer(model_tag):
    ## TODO check particle history:
    # Cannot check this right now, because reading the DPMJet event from ascii
    # randomly fails with a parsing error in HepMC3.

    ## TODO check vertices
    # Again, cannot check this right now, because it requires DPMJet events.
    # The situation is complicated by the fact that HepMC3 doesn't store vertices
    # of particles which do not have parents.

    generator = make_generator_instance(interaction_model_by_tag[model_tag])
    generator.init_generator(event_kinematics)

    test_file = "test_hepmc_writer_file_%s.dat" % model_tag

    event_data = []
    with HepMCWriter(test_file) as w:
        for event in generator.event_generator(event_kinematics, 3):
            event_data.append((event.pem_arr.T.copy(), event.p_ids.copy(), event.status.copy(), event.parents.T.copy() if event.parents is not None else 0))
            w.write(event)

    import pyhepmc_ng
    with pyhepmc_ng.open(test_file) as f:
        for ievent in range(3):
            event = f.read()
            assert event is not None
            assert event.event_number == ievent

            pem_ref, pid_ref, status_ref, parents_ref = event_data[ievent]

            particles = event.particles
            pem = np.empty((len(particles), 5))
            pid = np.empty(pem.shape[0], dtype=int)
            status = np.empty(pem.shape[0], dtype=int)
            for i, p in enumerate(particles):
                assert i + 1 == p.id
                pem[i] = p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e, p.generated_mass
                pid[i] = p.pid
                status[i] = p.status

            # parents = event.vertices
            # vt = np.empty((len(vertices), 4))
            # for i, v in enumerate(vertices):
            #     vt[i] = v.position.x, v.position.y, v.position.z, v.position.t

            assert_allclose(pem_ref, pem)
            assert_array_equal(pid_ref, pid)
            assert_array_equal(status_ref, status)
            # assert_array_equal(parents_ref, parents)
            # assert_allclose(vt_ref, vt)

    # delete test_file if test is successful
    import os
    os.unlink(test_file)
