"""
A commandline interface that mimics CRMC to ease transition from CRMC.

CRMC (Cosmic Ray Monte Carlo package) https://web.ikp.kit.edu/rulrich/crmc.html
"""

import argparse
import os
from . import models, __version__ as version
from .kinematics import CenterOfMass, FixedTarget, Momentum, _FromParticleName
from .util import AZ2pdg
from .constants import MeV, GeV
from .writer import Root
from pathlib import Path
from particle import Particle
from math import sqrt

MODELS = {
    # TODO add more numbers with our unique models
    # numbers must match definition in CRMC
    0: models.EposLHC,
    # 1: models.Epos199,
    2: models.QGSJet01d,
    # 3: models.Gheisha,
    4: models.Pythia6,
    # 5: models.Hijing138,
    6: models.Sibyll23d,
    7: models.QGSJetII04,
    8: models.Phojet193,
    11: models.QGSJetII03,
    # in CRMC 12 refers to different DPMJet versions
    12: models.DpmjetIII306,
    # models which are not in CRMC shoud use numbers upward of 20,
    # in case CRMC adds a few more
    20: models.Pythia8,
}
VALID_MODELS = ", ".join(f"{k}={v.label}" for (k, v) in MODELS.items())

FORMATS = (
    "hepmc",
    "hepmcgz",
    "root",
    # "lhe",
    # "lhegz",
)
VALID_FORMATS = f"{{{', '.join(FORMATS)}}}"


def extension(format):
    if format.endswith("gz"):
        return format[:-2] + ".gz"
    return format


def process_particle(x):
    try:
        x = int(x)
        if x > 10000:
            z = x // 10000
            a = (x % 10000) / 10
            return AZ2pdg(a, z)
        return x
    except ValueError:
        pass

    # handle any special names recognised by CRMC here
    # ...

    return _FromParticleName._get_pdg(x)


def particle_name(pid):
    return Particle.from_pdgid(pid).programmatic_name


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", action="store_true", help="print version")
    parser.add_argument(
        "-o", "--output", default="", help=f"output format {VALID_FORMATS}"
    )
    parser.add_argument(
        "-s",
        "--seed",
        type=int,
        default=0,
        help="seed between 1 and 2**31 (default is 0, which generates random seed)",
    )
    parser.add_argument(
        "-n",
        "--number",
        type=int,
        default=1,
        help="number of collisions (default is 1)",
    )
    parser.add_argument("-m", "--model", default=0, help=VALID_MODELS)
    parser.add_argument(
        "-p",
        "--projectile-momentum",
        type=float,
        default=0,
        help="projectile momentum in GeV/c",
    )
    parser.add_argument(
        "-P",
        "--target-momentum",
        type=float,
        default=0,
        help="target momentum in GeV/c",
    )
    parser.add_argument("-S", "--sqrts", type=float, default=0, help="sqrt(s) in GeV")
    parser.add_argument(
        "-i", "--projectile-id", default=2212, help="PDG ID or particle name"
    )
    parser.add_argument(
        "-I", "--target-id", default=2212, help="PDG ID or particle name"
    )
    parser.add_argument(
        "-f", "--out", help="Output file name (generated if none provided)"
    )

    args = parser.parse_args()

    if args.version:
        raise SystemExit(f"impy {version}")

    try:
        model_number = int(args.model)
        Model = MODELS[model_number]
    except KeyError:
        raise SystemExit(f"model {args.model} is invalid {VALID_MODELS}")
    except ValueError:
        matches = []
        x = args.model.lower()
        for M in MODELS.values():
            if x in M.label.lower():
                matches.append(M)
        if len(matches) == 1:
            Model = matches[0]
        else:
            raise SystemExit(f"model {args.model} is invalid {VALID_MODELS}")

    max_seed = int(1e9)
    if args.seed <= 0:
        args.seed = int.from_bytes(os.urandom(4), "little")
        args.seed %= max_seed
    else:
        if args.seed > max_seed:
            raise ValueError(f"{args.seed} is larger than maximum seed ({max_seed})")

    if args.number <= 0:
        raise SystemExit(f"invalid number of events {args.number}")

    pid1 = process_particle(args.projectile_id)
    pid2 = process_particle(args.target_id)

    p1 = args.projectile_momentum * GeV
    p2 = args.target_momentum * GeV
    sqrts = args.sqrts * GeV
    if sqrts < 0:
        raise SystemExit(f"sqrt(s) is negative {sqrts/GeV} GeV")
    if p1 != 0 or p2 != 0 and sqrts != 0:
        raise SystemExit("either set sqrts or momenta, not both")
    if p1 == 0 and p2 == 0 and sqrts == 0:
        raise SystemExit("no energy or momenta specified")

    if p1 != 0 and p2 != 0:
        # both momenta are non-zero, compute sqrt(s)
        m1 = Particle.from_pdgid(pid1).mass * MeV
        m2 = Particle.from_pdgid(pid2).mass * MeV
        e1 = sqrt(m1**2 + p1**2)
        e2 = sqrt(m2**2 + p2**2)
        # TODO use the numerically stable formula instead
        s = (e1 + e2) ** 2 - (p1 + p2) ** 2
        if s <= 0:
            raise SystemExit("s <= 0")
        sqrts = sqrt(s)

    if sqrts > 0:  # cms mode
        evt_kin = CenterOfMass(sqrts, pid1, pid2)
    elif p2 == 0:  # fixed target mode
        evt_kin = FixedTarget(Momentum(p1), pid1, pid2)

    if args.out == "-":
        args.output = "hepmc"
    else:
        if args.out:  # filename was provided
            args.out = Path(args.out)
            # try get format from filename extension
            format = "".join(x[1:] for x in args.out.suffixes[-2:])
            # check if both format and args.output are defined and are different
            if format and args.output and format != args.output:
                raise SystemExit(
                    f"File extension of {args.out} does not match {args.output}"
                )
            if not args.output:
                args.output = format or "hepmc"
        else:  # no filename provided
            if not args.output:
                args.output = "hepmc"
            # generate filename like CRMC
            ext = extension(args.output)
            p1 = particle_name(pid1)
            p2 = particle_name(pid2)
            mom = (sqrts or p1) / GeV
            fn = f"impy_{Model.pyname}_{args.seed}_{p1}_{p2}_{mom}.{ext}"
            odir = os.environ.get("CRMC_OUT", ".")
            args.out = Path(odir) / fn

        if not args.out.suffixes[-2:]:
            ext = "." + args.output
            if ext.endswith("gz"):
                ext = ext[:-2] + ".gz"
            args.out = Path(args.out).with_suffix(ext)

    if args.output not in FORMATS:
        raise SystemExit(f"unknown format {args.output} {VALID_FORMATS}")

    model = Model(evt_kin, seed=args.seed)

    if args.output.startswith("hepmc"):
        import pyhepmc

        writer = pyhepmc.open(args.out, "w")
    elif args.output == "root":
        writer = Root(args.out)

    with writer as w:
        for event in model(args.number):
            w.write(event)

    print(args.number, "events generated")


if __name__ == "__main__":
    main()
