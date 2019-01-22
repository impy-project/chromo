import sys
sys.path.append('../ParticleDataTool')
from impy.generators import epos
from impy.generators import sibyll
from impy.utils.units import TeV
from ParticleDataTool import PYTHIAParticleData
import numpy as np

pdg_table = PYTHIAParticleData()

def do_test(projectile, target, generator, version):
    nevent = 0
    for event in generator(version, projectile, target, 100):
        assert np.sum(event.px)
        assert np.sum(event.py)
        assert np.sum(event.pz)
        assert np.sum(event.en)
        assert np.sum(event.p_ids)
        

def test_epos_lhc():
    do_test(("p", 7*TeV), ("p", -7*TeV), epos, "LHC")

def test_sibyll_23c():
    do_test(("p", 7*TeV), ("p", -7*TeV), sibyll, "2.3c")