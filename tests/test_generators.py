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
from impy.util import get_all_models, name
import boost_histogram as bh
import gzip
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
from particle import literals as lp
import numpy as np
from scipy.stats import chi2

THIS_TEST = Path(__file__).stem
REFERENCE_PATH = Path(__file__).parent / "data" / THIS_TEST
REFERENCE_PATH.mkdir(exist_ok=True)


def run_model(Model, kin, number=1):
    try:
        gen = Model(kin, seed=1)
    except ValueError:
        return None

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

    values = []
    for _ in range(number):
        for event in gen(100):
            ev = event.final_state()
            h.fill(ev.eta, np.abs(ev.pid))
        if number == 1:
            return h
        else:
            values.append(h.values().copy())
            h.reset()
    return values


def compute_p_value(got, expected, expected_cov):
    delta = np.ravel(got) - np.ravel(expected)
    cov = np.reshape(expected_cov, delta.shape + delta.shape)
    for i in range(delta.size):
        cov[i, i] += 0.1  # prevent singularity
    inv_cov = np.linalg.inv(cov)
    v = np.einsum("i,ij,j", delta, inv_cov, delta)
    return 1 - chi2(delta.size).cdf(v)


def draw_comparison(fn, p_value, h, val_ref, cov_ref):
    fig, ax = plt.subplots(2, 3, figsize=(10, 8), constrained_layout=True)
    mname, projectile, target, frame = str(fn).split("_")
    plt.suptitle(f"{mname} {projectile} {target} {frame} pvalue={p_value}")
    xe = h.axes[0].edges
    cx = h.axes[0].centers
    for i, (axi, pdgid) in enumerate(zip(ax.flat, h.axes[1])):
        pname = name(pdgid)
        v = h.values()[:, i]
        vref = val_ref[:, i]
        cref = cov_ref[:, i, :, i]
        eref = np.diag(cref) ** 0.5
        vsum = np.sum(v)
        vrefsum = np.sum(vref)
        erefsum = np.einsum("i,ij,j", np.ones_like(v), cref, np.ones_like(v)) ** 0.5
        plt.sca(axi)
        plt.stairs(v, xe, color=f"C{i}", fill=True, alpha=0.5)
        plt.errorbar(cx, vref, eref, color=f"C{i}", marker="o")
        plt.title(f"{pname} {vsum:.0f} ({vrefsum:.0f} ± {erefsum:.0f})")
    fig_dir = Path() / THIS_TEST
    fig_dir.mkdir(exist_ok=True)
    plt.savefig(fig_dir / fn.with_suffix(".png"))
    plt.close(fig)


@pytest.mark.trylast
@pytest.mark.parametrize("target", ("p", "air"))
@pytest.mark.parametrize("projectile", ("gamma", "pi-", "p", "He"))
@pytest.mark.parametrize("frame", ("cms", "ft", "cms2ft", "ft2cms"))
@pytest.mark.parametrize("Model", get_all_models())
def test_generator(projectile, target, frame, Model):
    p1 = projectile
    p2 = target
    if p2 == "air":
        if Model is im.Pythia8:
            pytest.skip("Simulating nuclei in Pythia8 is very time-consuming")

        # cannot use Argon in SIBYLL, so make air from N, O only
        p2 = CompositeTarget((("N", 0.78), ("O", 0.22)))

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
        assert False  # we should never arrive here

    h = run_in_separate_process(run_model, Model, kin)
    if h is None:
        assert abs(kin.p1) not in Model.projectiles or abs(kin.p2) not in Model.targets
        return

    fn = Path(f"{Model.pyname}_{projectile}_{target}_{frame}")
    path_ref = REFERENCE_PATH / fn.with_suffix(".pkl.gz")
    if not path_ref.exists():
        print(f"{fn}: reference does not exist; generating... (but fail test)")
        # check plots to see whether reference makes any sense before committing it
        p_value = -1  # make sure test fails
        values = run_in_separate_process(run_model, Model, kin, 200, timeout=10000)
        val_ref = np.mean(values, axis=0)
        cov_ref = np.cov(np.transpose([np.ravel(x) for x in values]))
        cov_ref.shape = (*val_ref.shape, *val_ref.shape)
        with gzip.open(path_ref, "wb") as f:
            pickle.dump((val_ref, cov_ref), f)
    else:
        with gzip.open(path_ref) as f:
            val_ref, cov_ref = pickle.load(f)

        p_value = compute_p_value(h.values(), val_ref, cov_ref)

    draw_comparison(fn, p_value, h, val_ref, cov_ref)

    assert p_value >= 1e-3
