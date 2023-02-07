import subprocess as subp
from chromo import __version__ as version
from chromo import models as im
import re
from pathlib import Path
import pytest
from chromo.cli import MODELS
from particle import Particle
import pyhepmc
import uproot
import platform
import os


def format_matches_extension(p):
    ext = p.suffixes
    _, model, seed, pid1, pid2, en, *rest = p.stem.split("_")
    if ext[0] == ".hepmc":
        with pyhepmc.open(p) as f:
            event = f.read()
        assert event is not None
    elif ext[0] == ".root":
        with uproot.open(p) as f:
            tree = f["event"]
            assert tree.num_entries > 0


def run(
    *cmd,
    returncode=0,
    stdout=None,
    stderr=None,
    file=None,
    checks=(format_matches_extension,),
):
    r = subp.run(("chromo",) + cmd, capture_output=True)
    assert r.returncode == returncode, r.stderr.decode()
    match = None
    if stdout is not None:
        actual = r.stdout.decode()
        if stdout in actual:  # also check for literal match
            match = True
        else:
            match = re.search(stdout, actual)
        assert match, "\n" + actual
    if stderr is not None:
        actual = r.stderr.decode()
        if stderr in actual:  # also check for literal match
            match = True
        else:
            match = re.search(stderr, actual)
        assert match, f"\n  [expected] {stderr}\n  [  got   ] {actual}"
    if file:
        if match is not None:
            p = Path(file.format(*match.groups()))
        else:
            p = Path(file)
        assert p.exists()

        if checks:
            for check in checks:
                check(p)

        p.unlink()


# In the tests below, make sure that all output files
# have unique names, so that running tests in parallel does
# cause errors.


def test_no_args():
    run(
        returncode=1,
        stderr="Error: you need to specify particle momenta or the CMS energy",
    )


def test_version():
    run("-v", stdout=f"chromo {version}")


def test_minimum():
    run(
        "-S",
        "100",
        stdout="Seed[ \t]*([0-9]+)",
        file="chromo_eposlhc_{0}_2212_2212_100.hepmc",
    )


@pytest.mark.parametrize("seed", (0, -1))
def test_seed_1(seed):
    run(
        "-S",
        "100",
        "-s",
        f"{seed}",
        stdout="Seed[ \t]*([0-9]+)",
        file="chromo_eposlhc_{0}_2212_2212_100.hepmc",
    )


def test_seed_2():
    run(
        "-S",
        "100",
        "-s",
        "123",
        stdout="Seed[ \t]*123",
        file="chromo_eposlhc_123_2212_2212_100.hepmc",
    )


def test_number_1():
    run(
        "-S",
        "100",
        "-s",
        "1",
        "-n",
        "123",
        stdout="Collisions[ \t]*123",
        file="chromo_eposlhc_1_2212_2212_100.hepmc",
    )


def test_number_2():
    run("-S", "100", "-s", "1", "-n", "0", returncode=1)


@pytest.mark.parametrize(
    "spec,Model",
    tuple((str(k), v) for (k, v) in MODELS.items())
    + (
        ("eposlhc", im.EposLHC),
        ("sib23d", im.Sibyll23d),
        ("sib21", im.Sibyll21),
        ("sibyll-2.3", im.Sibyll23),
    ),
)
def test_model_1(spec, Model):
    run(
        "-S",
        "100",
        "-s",
        "2",
        "-m",
        spec,
        stdout=f"Model[ \t]*{Model.label}",
        file=f"chromo_{Model.pyname.lower()}_2_2212_2212_100.hepmc",
    )


@pytest.mark.skipif(
    platform.system() == "Windows", reason="Pythia-8 not available on Windows"
)
def test_model_2():
    run(
        "-m",
        "py8",
        returncode=1,
        stderr="Error: model=py8 is ambiguous, matches "
        "\\(Pythia-6.428, Pythia-8.[0-9]+\\)",
    )


def test_model_3():
    run(
        "-m",
        "sib",
        returncode=1,
        stderr=(
            "Error: model=sib is ambiguous, matches (SIBYLL-2.1, SIBYLL-2.3, "
            "SIBYLL-2.3c, SIBYLL-2.3d)"
        ),
    )


@pytest.mark.parametrize(
    "spec,expected",
    (
        ("pi+", 211),
        ("pi-", -211),
        ("p", 2212),
        ("O", 1000080160),
        ("211", 211),
        ("2212", 2212),
    ),
)
def test_projectile(spec, expected):
    p = Particle.from_pdgid(expected)
    name = p.name.replace("+", "\\+")
    run(
        "-S",
        "100",
        "-s",
        "3",
        "-i",
        spec,
        stdout=f"Projectile[ \t]*{name} \\({expected}\\)",
        file=f"chromo_eposlhc_3_{expected}_2212_100.hepmc",
    )


@pytest.mark.parametrize(
    "spec,expected",
    (
        ("p", 2212),
        ("C", 1000060120),
        ("O", 1000080160),
        ("2212", 2212),
    ),
)
def test_target(spec, expected):
    p = Particle.from_pdgid(expected)
    run(
        "-S",
        "100",
        "-s",
        "4",
        "-I",
        spec,
        stdout=f"Target[ \t]*{p.name} \\({expected}\\)",
        file=f"chromo_eposlhc_4_2212_{expected}_100.hepmc",
    )


def test_momentum_1():
    run(
        "-s",
        "5",
        "-p 1000",
        "-P -1000",
        stdout="sqrt\\(s\\)[ \t]*2000 GeV",
        file="chromo_eposlhc_5_2212_2212_2000.hepmc",
    )


def test_momentum_2():
    run(
        "-s",
        "6",
        "-p 1000",
        "-P 0",
        stdout="""\
[ |]*Projectile momentum[ \t]*1000 GeV/c *|
[ |]*Target momentum[ \t]*0 GeV/c *|
[ |]*sqrt\\(s\\)[ \t]*43.3394 GeV\
""",
        file="chromo_eposlhc_6_2212_2212_43.hepmc",
    )


def test_format_1():
    run(
        "-s",
        "7",
        "-S",
        "100",
        stdout="Format[ \t]*hepmc",
        file="chromo_eposlhc_7_2212_2212_100.hepmc",
    )


@pytest.mark.parametrize("format", ("hepmc", "hepmcgz", "root"))
@pytest.mark.parametrize("model", ("EPOS-LHC", "SIBYLL-2.1", "Pythia-6.4"))
def test_format_2(format, model):
    ext = format
    if ext.endswith("gz"):
        ext = ext[:-2] + ".gz"

    pyname = {"EPOS-LHC": "eposlhc", "SIBYLL-2.1": "sibyll21", "Pythia-6.4": "pythia6"}[
        model
    ]

    run(
        "-s",
        "8",
        "-S",
        "100",
        "-o",
        format,
        "-m",
        model,
        stdout=f"Format[ \t]*{format}",
        file=f"chromo_{pyname}_8_2212_2212_100.{ext}",
    )


@pytest.mark.skipif(
    "CI" in os.environ,
    reason="skip on CI, because cibuildwheel has problems with "
    "graphviz installation",
)
def test_format_3():
    if platform.system() == "Windows":
        pytest.xfail(
            "Test aborts on Windows with this message: "
            "UnicodeEncodeError: 'charmap' codec can't encode character '\u0394'"
            " in  position 20049: character maps to <undefined>"
        )
    run(
        "-s",
        "9",
        "-S",
        "100",
        "-o",
        "svg",
        stdout="Format[ \t]*svg",
        file="chromo_eposlhc_9_2212_2212_100_000.svg",
    )
