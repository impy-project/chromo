from impy import util
import numpy as np
from numpy.testing import assert_equal
from particle import literals as lp


def test_select_parents():
    mask = np.array([False, True, False])
    assert util.select_parents(mask, None) is None

    parents = np.array([[0, 0], [1, 2], [2, 0]])
    par = util.select_parents(mask, parents)
    assert_equal(par, [[0, 0]])

    mask = np.array([False, True, True])
    parents = np.array([[0, 0], [1, 2], [2, 0]])
    par = util.select_parents(mask, parents)
    assert_equal(par, [[0, 0], [1, 0]])


def test_tolerant_string_match():
    assert util.tolerant_string_match("123", "abc1 s2.d34")
    assert not util.tolerant_string_match("123", "321")
    assert not util.tolerant_string_match("123", "3123")
    assert not util.tolerant_string_match("123", "124")


def test_Nuclei():
    photon = lp.photon.pdgid
    proton = lp.proton.pdgid
    carbon = util.name2pdg("C")
    lead = util.name2pdg("Pb")

    a = util.Nuclei()
    assert photon not in a
    assert proton in a
    assert carbon in a
    assert lead in a
    assert repr(a) == "Nuclei(a_min=1, a_max=1000, z_min=0, z_max=1000)"

    b = {photon} | a
    assert photon in b
    assert proton in b
    assert carbon in b
    assert lead in b
    assert repr(b) == "Nuclei(a_min=1, a_max=1000, z_min=0, z_max=1000) | {<PDGID: 22>}"

    c = a | {photon}
    assert photon in c
    assert proton in c
    assert carbon in c
    assert lead in c
    assert repr(c) == "Nuclei(a_min=1, a_max=1000, z_min=0, z_max=1000) | {<PDGID: 22>}"

    d = util.Nuclei(a_max=20)
    assert photon not in d
    assert proton in d
    assert carbon in d
    assert lead not in d
    assert repr(d) == "Nuclei(a_min=1, a_max=20, z_min=0, z_max=1000)"
