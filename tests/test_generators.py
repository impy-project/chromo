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
from impy.util import get_all_models, pdg2name
import impy
import boost_histogram as bh
import gzip
import pickle
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
from particle import literals as lp
import numpy as np
from scipy.stats import chi2
import sys

matplotlib.use("svg")  # need non-interactive backend for CI Windows

THIS_TEST = Path(__file__).stem
REFERENCE_PATH = Path(__file__).parent / "data" / THIS_TEST
REFERENCE_PATH.mkdir(exist_ok=True)
FIG_PATH = Path("fig")
FIG_PATH.mkdir(parents=True, exist_ok=True)


def run_model(Model, kin, number=1):
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

    def fill(gen, kin, h):
        for event in gen(kin, 100):
            ev = event.final_state()
            h.fill(ev.eta, np.abs(ev.pid))

    try:
        if number == 1:
            gen = Model(seed=1)
            fill(gen, kin, h)
            return h
        else:
            values = []
            gen = Model(seed=1)
            for i in range(number):
                sys.stderr.write(f"iteration {i}")
                sys.stderr.flush()
                fill(gen, kin, h)
                values.append(h.values().copy())
                h.reset()
            return values
    except ValueError:
        return None


def compute_p_value(got, expected, cov):
    delta = got - expected
    var = np.diag(cov) + 0.1  # prevent singularity
    # we ignore correlations, since including them
    # leads to implausible results
    v = np.sum(delta**2 / var)
    ndof = np.sum(delta != 0)
    return 1 - chi2(ndof).cdf(v)


def draw_comparison(fn, p_value, axes, values, val_ref, cov_ref):
    import re

    m = re.search(r"(.+)_(.+)_(.+)_(.+)", str(fn))
    mname = m.group(0)
    projectile = m.group(1)
    target = m.group(2)
    frame = m.group(3)
    fig = plt.figure(figsize=(10, 8))
    plt.suptitle(f"{mname} {projectile} {target} {frame} pvalue={p_value:.3g}")
    xe = axes[0].edges
    cx = axes[0].centers
    left = 0.1
    bottom = 0.07
    width = 0.25
    height1 = 0.25
    height2 = 0.1
    hspace = 0.05
    wspace = 0.05
    shape = tuple(len(x) for x in axes)
    val = np.reshape(values, shape)
    for i, pdgid in enumerate(axes[1]):
        pname = pdg2name(pdgid)
        v = val[:, i]
        vref = np.reshape(val_ref, shape)[:, i]
        cref = np.reshape(cov_ref, shape * 2)[:, i, :, i]
        eref = np.diag(np.maximum(cref, 0)) ** 0.5
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
        plt.xlim(xe[0], xe[-1])
        plt.tick_params(bottom=False, labelbottom=False)
        ax = fig.add_axes(
            (
                left + (width + wspace) * icol,
                bottom + (1 - irow) * (height1 + height2 + 2 * hspace),
                width,
                height2,
            )
        )
        plt.sca(ax)
        d = (v - vref) / (eref + 0.1)
        plt.plot(cx, d, color="k")
        plt.ylim(-5, 5)
        plt.xlim(xe[0], xe[-1])
    plt.savefig(FIG_PATH / fn.with_suffix(".svg"))
    plt.close(fig)


@pytest.mark.trylast
@pytest.mark.parametrize("frame", ("cms", "ft", "cms2ft", "ft2cms"))
@pytest.mark.parametrize("target", ("p", "air"))
@pytest.mark.parametrize("projectile", ("gamma", "pi-", "p", "He"))
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
        kin = EventKinematics(p1, p2, ecm=100 * GeV, frame=EventFrame.FIXED_TARGET)
    elif frame == "ft2cms":
        kin = EventKinematics(p1, p2, elab=100 * GeV, frame=EventFrame.CENTER_OF_MASS)
    else:
        assert False  # we should never arrive here

    h = run_model(Model, kin)
    if h is None:
        assert abs(kin.p1) not in Model.projectiles or abs(kin.p2) not in Model.targets
        return

    fn = Path(f"{Model.pyname}_{projectile}_{target}_{frame}")
    path_ref = REFERENCE_PATH / fn.with_suffix(".pkl.gz")
    if not path_ref.exists():
        print(f"{fn}: reference does not exist; generating...")
        # check plots to see whether reference makes any sense before committing it
        values = run_model(Model, kin, 50)
        val_ref = np.reshape(np.mean(values, axis=0), -1)
        cov_ref = np.cov(np.transpose([np.reshape(x, -1) for x in values]))
        assert np.sum(np.isnan(val_ref)) == 0
        assert np.sum(np.isnan(cov_ref)) == 0
        with gzip.open(path_ref, "wb") as f:
            pickle.dump((val_ref, cov_ref), f)
        values = np.reshape(h.values(), -1)
        draw_comparison(fn, -1, h.axes, values, val_ref, cov_ref)
        pytest.xfail(reason="reference does not exist; generated new one, check it")

    with gzip.open(path_ref) as f:
        val_ref, cov_ref = pickle.load(f)

    # histogram sometimes contains extreme spikes
    # not reflected in cov, this murky formula
    # protects against these false negatives
    values = np.reshape(h.values(), -1)
    for i, v in enumerate(values):
        cov_ref[i, i] = 0.5 * (cov_ref[i, i] + v)

    p_value = compute_p_value(values, val_ref, cov_ref)

    threshold = 1e-6

    if not (p_value >= threshold) or impy.debug_level > 0:
        draw_comparison(fn, p_value, h.axes, values, val_ref, cov_ref)

    assert p_value >= threshold
