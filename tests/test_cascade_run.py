import sys
sys.path.append('../ParticleDataTool')
from impy import common, sibyll
from impy.utils import EnergySpectrumBoost
import numpy as np

def test_sibyll():
    sibyll_run = sibyll.SibyllCascadeRun(
        lib_str = 'sib23c',
        label = 'SIBYLL 2.3c',
        decay_mode = 2,
        n_events=1000,
        fill_subset=True,
        evt_class=sibyll.SibyllCascadeEvent)

    evkin = common.EventKinematics(plab=1e8, p1pdg=2212, nuc2_prop=(14,7))

    spectrum_hists = []
    pdt = sibyll_run.ptab
    for pdgid in [211, -211, 321, -321, 2212, 2112]:
        sibyll_run.spectrum_hists.append(
            EnergySpectrumBoost(lorentz_factors=(evkin.gamma_cm, evkin.betagamma_cm),
                            part_id = pdt.pdg2modid[pdgid],
                            pdg_id = pdgid,
                            E_lab=evkin.elab,
                            title = pdt.pdg2modname[pdgid],
                            grid_var = 'x',
                            x_bins=np.logspace(-5,1,100))
        )

    sibyll_run.start('p', evkin.elab, evkin.ecm, Atarget=1)

    for h in sibyll_run.spectrum_hists:
        assert sum(h.get()) > 0