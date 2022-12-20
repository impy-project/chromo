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


def compute_p_value(got, expected, cov):
    delta = np.reshape(got, -1) - expected
    for i in range(delta.size):
        cov[i, i] += 0.1  # prevent singularity
    inv_cov = np.linalg.inv(cov)
    v = np.einsum("i,ij,j", delta, inv_cov, delta)
    ndof = np.sum(delta != 0)
    return 1 - chi2(ndof).cdf(v)


def draw_comparison(fn, p_value, h, val_ref, cov_ref):
    mname, projectile, target, frame = str(fn).split("_")
    fig = plt.figure(figsize=(10, 8))
    plt.suptitle(f"{mname} {projectile} {target} {frame} pvalue={p_value:.3g}")
    xe = h.axes[0].edges
    cx = h.axes[0].centers
    left = 0.1
    bottom = 0.07
    width = 0.25
    height1 = 0.25
    height2 = 0.1
    hspace = 0.05
    wspace = 0.05
    val = h.values()
    for i, pdgid in enumerate(h.axes[1]):
        pname = name(pdgid)
        v = h.values()[:, i]
        vref = np.reshape(val_ref, val.shape)[:, i]
        cref = np.reshape(cov_ref, val.shape * 2)[:, i, :, i]
        eref = np.diag(cref) ** 0.5
        vsum = np.sum(v)
        vrefsum = np.sum(vref)
        erefsum = np.sum(cref) ** 0.5
        icol = i % 3
        irow = i // 3
        ax = fig.add_axes(
            (
                left + (width + wspace) * icol,
                bottom
                + (1 - irow) * (height1 + height2 + 2 * hspace)
                + height2
                + 0.5 * hspace,
                width,
                height1,
            )
        )
        plt.sca(ax)
        plt.stairs(v, xe, color=f"C{i}", fill=True, alpha=0.5)
        plt.errorbar(cx, vref, eref, color=f"C{i}", marker="o")
        plt.title(f"{pname} {vsum:.0f} ({vrefsum:.0f} ± {erefsum:.0f})")
        plt.tick_params(bottom=False, labelbottom=False)
        ax = fig.add_axes(
            (
                left + (width + wspace) * icol,
                bottom + (1 - irow) * (height1 + height2 + 2 * hspace),
                width,
                height2,
            )
        )
        plt.xlim(xe[0], xe[-1])
        plt.sca(ax)
        d = (v - vref) / (eref + 0.1)
        plt.plot(cx, d, color="k")
        plt.ylim(-5, 5)
        plt.xlim(xe[0], xe[-1])
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
        values = run_in_separate_process(run_model, Model, kin, 50, timeout=10000)
        val_ref = np.reshape(np.mean(values, axis=0), -1)
        cov_ref = np.cov(np.transpose([np.reshape(x, -1) for x in values]))
        with gzip.open(path_ref, "wb") as f:
            pickle.dump((val_ref, cov_ref), f)
    else:
        with gzip.open(path_ref) as f:
            val_ref, cov_ref = pickle.load(f)

        p_value = compute_p_value(h.values(), val_ref, cov_ref)

    draw_comparison(fn, p_value, h, val_ref, cov_ref)

    assert p_value >= 1e-3
