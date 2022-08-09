import numpy as np

from impy.constants import GeV
from impy.kinematics import EventKinematics

# from impy.common import impy_config
from impy.models import UrQMD34


def test_urqmd():

    # AF: This is what the user interaction has to yield.
    # It is the typical expected configuration that one
    # wants to run (read pp-mode at energies not exceeding
    # 7 TeV). If you want cosmic ray energies, this should
    # be rather p-N at 10 EeV and lab frame (not yet defined).

    event_kinematics = EventKinematics(ecm=200 * GeV, p1pdg=2112, nuc2_prop=(14, 7))

    # impy_config["user_frame"] = 'laboratory'
    # impy_config['tau_stable'] = 1.
    # impy_config['pi0_stable'] = False
    generator = UrQMD34(event_kinematics)
    # import IPython
    # IPython.embed()

    make_decay_list = [211, 321, 111, 2112, 310, 130, 13, -13, 3122, -3212]
    for pid in make_decay_list:
        generator.set_stable(pid, stable=False)
    generator.lib.stables.stabvec *= 0
    generator.lib.stables.nstable = 0

    egrid = np.geomspace(1e-1 * GeV, event_kinematics.ecm * 1.2, 30)
    pions = 0
    protons = 0

    for event in generator(10):
        protons += np.histogram(event.en[event.p_ids == 2212], bins=egrid)[0]
        pions += np.histogram(event.en[np.abs(event.p_ids) == 211], bins=egrid)[0]
        # generator.lib.pydat3.mdcy[102 - 1, 0] = 1
        # event.filter_final_state()
        # print 'px', event.px
        # print 'py', event.py
        # print 'pz', event.pz
        # print 'en', event.en
        # print 'p_ids', event.p_ids
        # print 'impact param', event.impact_parameter
        # import IPython
        # IPython.embed()
        # print event.impact_parameter, event.n_wounded_A, event.n_wounded_B#,
        # event.n_NN_interactions

    # TODO better tests
    assert np.sum(pions) > 0
    assert np.sum(protons) > 0

    # AF: Maybe it would be better, or a good alternative
    # to make a Particle class and the MCEvent provides an iterator
    # to it's contents as objects of this Particle class. This would
    # be Pythia 8 style but probably slow in PYTHON.

    # for event in ....:
    #     for particle in event:
    #         print particle.p_id, particle.en...
