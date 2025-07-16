import numpy as np
from numpy.testing import assert_equal
from particle import literals as lp
from pytest import approx

from chromo import util


def test_select_mothers():
    mask = np.array([False, True, False])
    assert util.select_mothers(mask, None) is None

    mothers = np.array([[-1, -1], [0, 1], [1, -1]])
    par = util.select_mothers(mask, mothers)
    assert_equal(par, [[-1, -1]])

    mask = np.array([False, True, True])
    mothers = np.array([[-1, -1], [0, 1], [1, -1]])
    par = util.select_mothers(mask, mothers)
    assert_equal(par, [[-1, -1], [0, -1]])


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


def test_ecm_elab_conversion():
    elab = 2.1
    m1 = 1.1
    m2 = 0.5
    ecm = util.elab2ecm(elab, m1, m2)
    elab2 = util.ecm2elab(ecm, m1, m2)
    assert elab2 == approx(elab)


def test_momentum_energy_conversion():
    p = 2.1
    m = 1.1
    en = util.momentum2energy(p, m)
    p2 = util.energy2momentum(en, m)
    assert p == approx(p2)


def test_is_real_nucleus():
    proton = 2212
    neutron = 2112

    assert not util.is_real_nucleus(proton)
    assert not util.is_real_nucleus(neutron)

    proton2 = util.AZ2pdg(1, 1)
    neutron2 = util.AZ2pdg(1, 0)

    assert not util.is_real_nucleus(proton2)
    assert not util.is_real_nucleus(neutron2)

    deuterium = util.AZ2pdg(2, 1)

    assert util.is_real_nucleus(deuterium)

    mix = util.CompositeTarget([("p", 0.5), ("He", 0.5)])

    assert util.is_real_nucleus(mix)
