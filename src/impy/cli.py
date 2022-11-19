"""
A commandline interface that mimics CRMC to ease transition from CRMC.

CRMC (Cosmic Ray Monte Carlo package) https://web.ikp.kit.edu/rulrich/crmc.html
"""

import argparse
import os
from . import models, __version__ as version
from .kinematics import CenterOfMass, FixedTarget, Momentum, _FromParticleName
from .util import AZ2pdg, tolerant_string_match, get_all_models
from .constants import MeV, GeV
from .writer import Root
from pathlib import Path
from particle import Particle
from math import sqrt

# Only add numbers here of models that to match CRMC.
# Impy models which are not in CRMC should not be added here.
MODELS = {
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
}
VALID_MODELS = ", ".join(f"{k}={v.label}" for (k, v) in MODELS.items())
VALID_MODELS += ", ".join(
    {v.label for v in get_all_models() if v not in MODELS.values()}
)

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


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", action="store_true", help="print version")
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
    parser.add_argument(
        "-m",
        "--model",
        default=0,
        help=f"select model via number or name with tolerant matching {VALID_MODELS}",
    )
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
        "-i",
        "--projectile-id",
        default=2212,
        help="PDG ID or particle name, default is proton",
    )
    parser.add_argument(
        "-I",
        "--target-id",
        default=2212,
        help="PDG ID or particle name, default is proton",
    )
    parser.add_argument(
        "-o", "--output", default="", help=f"output format {VALID_FORMATS}"
    )
    parser.add_argument(
        "-f", "--out", help="Output file name (generated if none provided)"
    )

    args = parser.parse_args()

    if args.version:
        raise SystemExit

    max_seed = int(1e9)  # EPOS requirement
    if args.seed <= 0:
        args.seed = int.from_bytes(os.urandom(4), "little")
        args.seed %= max_seed
    else:
        if args.seed > max_seed:
            raise SystemExit(
                f"Error: {args.seed} is larger " f"than maximum seed ({max_seed})"
            )

    if args.number <= 0:
        raise SystemExit("Error number must be positive")

    try:
        model_number = int(args.model)
        Model = MODELS[model_number]
    except KeyError:
        raise SystemExit(f"Error: model {args.model} is invalid {VALID_MODELS}")
    except ValueError:
        matches = []
        x = args.model.lower()
        for M in get_all_models():
            if tolerant_string_match(x, M.label.lower()):
                matches.append(M)
        if len(matches) == 1:
            Model = matches[0]
        elif len(matches) == 0:
            raise SystemExit(f"Error: model {args.model} is invalid {VALID_MODELS}")
        else:
            raise SystemExit(
                f"Error: model {args.model} is ambiguous, "
                f"matches {', '.join(v.label for v in matches)}"
            )

    args.projectile_id = process_particle(args.projectile_id)
    args.target_id = process_particle(args.target_id)

    args.projectile_momentum *= GeV
    args.target_momentum *= GeV
    args.sqrts *= GeV
    pr = args.projectile_momentum
    ta = args.target_momentum
    if args.sqrts < 0:
        raise SystemExit(f"Error: sqrt(s) is negative {args.sqrts/GeV} GeV")
    if (pr != 0 or ta != 0) and args.sqrts != 0:
        raise SystemExit("Error: either set sqrts or momenta, but not both")
    if pr == 0 and ta == 0 and args.sqrts == 0:
        raise SystemExit("Error: no energy (-S) or momenta specified (-p, -P)")

    if pr != 0 or ta != 0:
        # compute sqrt(s)
        m1 = Particle.from_pdgid(args.projectile_id).mass * MeV
        m2 = Particle.from_pdgid(args.target_id).mass * MeV
        e1 = sqrt(m1**2 + pr**2)
        e2 = sqrt(m2**2 + ta**2)
        # TODO use the numerically stable formula instead
        s = (e1 + e2) ** 2 - (pr + ta) ** 2
        if s <= 0:
            raise SystemExit("Error: s <= 0")
        args.sqrts = sqrt(s)

    if args.out == "-":
        args.output = "hepmc"
        raise SystemExit("Error: output to stdout not yet supported")
    else:
        if args.out:  # filename was provided
            args.out = Path(args.out)
            # try get format from filename extension
            format = "".join(x[1:] for x in args.out.suffixes[-2:])
            # check if both format and args.output are defined and are different
            if format and args.output and format != args.output:
                raise SystemExit(
                    f"Error: File extension of {args.out} does not match {args.output}"
                )
            if not args.output:
                args.output = format or "hepmc"
        else:  # no filename provided
            if not args.output:
                args.output = "hepmc"
            # generate filename like CRMC
            ext = extension(args.output)
            pid1 = args.projectile_id
            pid2 = args.target_id
            en = f"{args.sqrts / GeV:.0f}"
            fn = f"impy_{Model.pyname.lower()}_{args.seed}_{pid1}_{pid2}_{en}.{ext}"
            odir = os.environ.get("CRMC_OUT", ".")
            args.out = Path(odir) / fn

        if not args.out.suffixes[-2:]:
            ext = "." + args.output
            if ext.endswith("gz"):
                ext = ext[:-2] + ".gz"
            args.out = Path(args.out).with_suffix(ext)

    if args.output not in FORMATS:
        raise SystemExit(f"Error: unknown format {args.output} {VALID_FORMATS}")

    args.model = Model
    return args


def main():
    print(f"impy {version}")

    args = parse_arguments()

    print(f"  Model: {args.model.label}")
    print(f"  Collisions: {args.number}")
    p = Particle.from_pdgid(args.projectile_id)
    print(f"  Projectile: {p.name} ({args.projectile_id})")
    p = Particle.from_pdgid(args.target_id)
    print(f"  Target: {p.name} ({args.target_id})")
    print(f"  Seed: {args.seed}")
    pr = args.projectile_momentum
    ta = args.target_momentum
    if pr != 0 or ta != 0:
        print(f"  Projectile momentum: {pr:g} GeV/c")
        print(f"  Target momentum: {ta:g} GeV/c")
    print(f"  sqrt(s): {args.sqrts:g} GeV")
    print(f"  Format: {args.output}")

    if pr > 0 and ta == 0:  # fixed target mode
        evt_kin = FixedTarget(Momentum(pr), args.projectile_id, args.target_id)
    else:  # cms mode
        evt_kin = CenterOfMass(args.sqrts, args.projectile_id, args.target_id)

    if args.output.startswith("hepmc"):
        import pyhepmc

        writer = pyhepmc.open(args.out, "w")
    elif args.output == "root":
        writer = Root(args.out)

    model = args.model(evt_kin)

    with writer as w:
        for event in model(args.number):
            w.write(event)


if __name__ == "__main__":
    main()
