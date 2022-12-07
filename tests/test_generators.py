from impy.constants import GeV
from impy.kinematics import (
    CenterOfMass,
    EventKinematics,
    FixedTarget,
    TotalEnergy,
    EventFrame,
    CompositeTarget,
)
import impy.models as im
import pytest
from .util import run_in_separate_process
from impy.util import get_all_models
import boost_histogram as bh
import gzip
import pickle
import os
from pathlib import Path
from numpy.testing import assert_allclose
from particle import Particle, literals as lp
import numpy as np

# generate list of all models in impy.models
models = get_all_models()
ref_dir = Path(__file__).parent / "data"


def run_model(Model, evt_kin):
    gen = Model(evt_kin, seed=1)

    h = bh.Histogram(
        bh.axis.Regular(21, -10, 10),
        bh.axis.IntCategory(
            [
                p.pdgid
                for p in (
                    lp.pi_plus,
                    lp.K_plus,
                    lp.K_S_0,
                    lp.K_L_0,
                    lp.proton,
                    lp.neutron,
                )
            ]
        ),
    )
    for event in gen(100):
        ev = event.final_state()
        h.fill(ev.eta, np.abs(ev.pid))

    return h


@pytest.mark.parametrize("target", ("p", "air"))
@pytest.mark.parametrize("frame", ("cms", "ft", "cms2ft", "ft2cms"))
@pytest.mark.parametrize("Model", models)
def test_generator(target, frame, Model):
    p1 = lp.pi_minus.pdgid
    if Model is im.Sophia20:
        # Sophia can only do γp, γn
        p1 = lp.gamma.pdgid
    elif Model in [im.Phojet112, im.UrQMD34]:
        # The old phojet needs more tweaking for pion-proton (is not related to test)
        p1 = lp.proton.pdgid

    if target == "air":
        if Model is im.Pythia8:
            pytest.skip("Simulating nuclei in Pythia8 is very time-consuming")

        # cannot use Argon in SIBYLL, so make air from N, O only
        p2 = CompositeTarget((("N", 0.78), ("O", 0.22)))
        for c in p2.components:
            if c not in Model.targets:
                pytest.skip(f"{Model.pyname} does not support support nuclei")
    else:
        p2 = target

    if frame == "cms":
        kin = CenterOfMass(100 * GeV, p1, p2)
    elif frame == "ft":
        kin = FixedTarget(TotalEnergy(100 * GeV), p1, p2)
    elif frame == "cms2ft":
        kin = EventKinematics(
            ecm=100 * GeV, particle1=p1, particle2=p2, frame=EventFrame.FIXED_TARGET
        )
    elif frame == "ft2cms":
        kin = EventKinematics(
            elab=100 * GeV, particle1=p1, particle2=p2, frame=EventFrame.CENTER_OF_MASS
        )
    else:
        assert False

    h = run_in_separate_process(run_model, Model, kin)

    path = Path(f"test_generators_{Model.pyname}_{target}_{frame}.pkl.gz")
    try:
        path_ref = ref_dir / path
        with gzip.open(path_ref) as f:
            href = pickle.load(f)
        assert_allclose(h.values(), href.values(), atol=1, rtol=0.05)
    except (AssertionError, FileNotFoundError) as exc:
        if "CI" not in os.environ:
            # when run locally, generate plots for visual inspection if test fails
            try:
                import matplotlib.pyplot as plt

                _, ax = plt.subplots(1, 2, figsize=(10, 5), sharex=True, sharey=True)
                plt.suptitle(f"{Model.pyname} {target} {frame}")
                for i, pdgid in enumerate(h.axes[1]):
                    p = Particle.from_pdgid(pdgid)
                    v = h.values()[:, i]
                    xe = h.axes[0].edges
                    ax[0].stairs(v, xe, label=f"{p.name} {np.sum(v):.0f}")
                    if not isinstance(exc, FileNotFoundError):
                        vref = href.values()[:, i]
                        ax[1].stairs(vref, xe, label=f"{p.name}{np.sum(vref):.0f}")
                    for axi in ax:
                        axi.legend(
                            loc="upper center",
                            title="c.c. included",
                            ncol=3,
                            frameon=False,
                        )
                plt.semilogy()
                plt.savefig(path.with_suffix(".png"))

            except ModuleNotFoundError:
                print("I want to show some plots, please install matplotlib")

            with gzip.open(path, "wb") as f:
                pickle.dump(h, f)

        raise
