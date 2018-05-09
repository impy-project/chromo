import sys
import os
import numpy as np

sys.path.append(os.path.dirname(__file__) + "/..")

from impy.constants import *
from impy.kinematics import EventKinematics
from impy.models.sibyll import SIBYLLRun
from impy.models.dpmjetIII import DpmjetIIIRun
from impy.common import impy_config, pdata
from impy.writer import HepMCWriter

event_kinematics = EventKinematics(
    ecm=7 * TeV,
    p1pdg=2212,
    nuc2_prop=(16,8))

impy_config["user_frame"] = 'laboratory'


def test_hepevt_access():
    import dpmjet306
    generator = DpmjetIIIRun(dpmjet306)
    generator.init_generator(event_kinematics)

    import pyhepmc
    for event in generator.event_generator(event_kinematics, 1):
        # print 'pz', event.pz
        # print 'p_ids', event.p_ids
        hepevt = event.lib.dtevt1
        import IPython
        IPython.embed()
        # print typ
        # print hepevt.nhkk
        # print hepevt.nevhkk
        # address = event.lib.dtevt1.__array_interface__['data'][0]
        # pyhepmc.print_hepevt(address)


def test_hepmc_writer():
    import dpmjet306
    generator = DpmjetIIIRun(dpmjet306)
    generator.init_generator(event_kinematics)

    with HepMCWriter("testfile.dat") as w:
        for event in generator.event_generator(event_kinematics, 1):
            w.write(event)


if __name__ == '__main__':
    # test_hepevt_access()
    test_hepmc_writer()