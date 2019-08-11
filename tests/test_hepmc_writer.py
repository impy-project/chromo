from __future__ import print_function

import sys
import os
import numpy as np

root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))
sys.path.append(root_dir)

from impy.definitions import *
from impy.constants import *
from impy.kinematics import EventKinematics
from impy.common import impy_config, pdata
from impy.util import info

event_kinematics = EventKinematics(ecm=7000 * GeV,
                                   p1pdg=2212,
                                   p2pdg=2212)

impy_config["user_frame"] = 'laboratory'

# def test_hepevt_access():
#     generator = make_generator_instance(interaction_model_by_tag['SIBYLL23C'])
#     generator.init_generator(event_kinematics)
# 
#     for event in generator.event_generator(event_kinematics, 1):
#         # print 'pz', event.pz
#         # print 'p_ids', event.p_ids
#         hepevt = event.lib.dtevt1
#         import IPython
#         IPython.embed()
#         # print typ
#         # print hepevt.nhkk
#         # print hepevt.nevhkk
#         # address = event.lib.dtevt1.__array_interface__['data'][0]
#         # pyhepmc.print_hepevt(address)


def test_hepmc_writer():
    from impy.writer import HepMCWriter

    generator = make_generator_instance(interaction_model_by_tag['SIBYLL23C'])
    generator.init_generator(event_kinematics)

    test_file = "test_hepmc_writer_file.dat"
    with HepMCWriter(test_file) as w:
        for event in generator.event_generator(event_kinematics, 1):
            w.write(event)

    with open(test_file) as f:
        content = f.read()
        assert content.startswith("HepMC::Version")
        assert content.endswith("END_EVENT_LISTING\n\n")

    # delete test_file if test is successful
    import os
    os.unlink(test_file)


if __name__ == '__main__':
    # test_hepevt_access()
    test_hepmc_writer()