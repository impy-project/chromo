"""
A commandline interface that mimics CRMC to ease transition from CRMC.

CRMC (Cosmic Ray Monte Carlo package) https://web.ikp.kit.edu/rulrich/crmc.html
"""

import argparse
import os
from chromo import models, __version__ as version
from chromo.kinematics import CenterOfMass, FixedTarget, Momentum
from chromo.util import AZ2pdg, tolerant_string_match, get_all_models, name2pdg
from chromo.constants import MeV, GeV
from chromo import writer
from pathlib import Path
from particle import Particle
from math import sqrt
from rich.progress import (
    Progress,
    ProgressColumn,
    Task,
    BarColumn,
    MofNCompleteColumn,
    TimeRemainingColumn,
    TaskProgressColumn,
)
from rich.text import Text


class SpeedColumn(ProgressColumn):
    """Renders generation speed."""

    def render(self, task: "Task") -> Text:
        """Show generation speed."""
        speed = task.finished_speed or task.speed
        if speed is None:
            return Text("?", style="progress.data.speed")
        return Text(f"{speed:.0f}/s", style="progress.data.speed")


# Only add numbers here for backward-compatibility with CRMC.
# Chromo-exclusive models do not get a number and should not be added here.
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
    # 12 refers to different DPMJet versions in CRMC
    12: models.DpmjetIII307,
}

VALID_MODELS = []
for M in get_all_models():
    for k, v in MODELS.items():
        if v is M:
            VALID_MODELS.append(f"{M.label} [{k}]")
            break
    else:
        VALID_MODELS.append(M.label)
VALID_MODELS = ", ".join(sorted(VALID_MODELS))

FORMATS = {
    "hepmc": writer.Hepmc,
    "hepmc:gz": writer.Hepmc,
    "root": writer.Root,
    "root:vertex": lambda *args: writer.Root(*args, write_vertices=True),
    "svg": writer.Svg,
    "null": writer.Null,
    # "lhe",
    # "lhegz",
}
VALID_FORMATS = f"{', '.join(FORMATS)}"


def extension(format):
    extension, *options = format.split(":")
    if "gz" in options:
        return f"{extension}.gz"
    return extension


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

    try:
        return name2pdg(x)
    except KeyError:
        raise SystemExit(f"particle name {x} not recognized")


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", action="store_true", help="print version")
    parser.add_argument(
        "-s",
        "--seed",
        type=int,
        default=0,
        help="seed between 1 and 1e9 (default is 0, which generates random seed)",
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
        help=(
            f"select model via number or via tolerant string match (case-insensitive), "
            f"allowed values: {VALID_MODELS}"
        ),
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
        default="p",
        help="PDG ID or particle name (default is proton)",
    )
    parser.add_argument(
        "-I",
        "--target-id",
        default="p",
        help="PDG ID or particle name (default is proton)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="",
        help=f"output format: {VALID_FORMATS} (default is hepmc)",
    )
    parser.add_argument(
        "-f",
        "--out",
        help="output file name (generated if none provided); "
        "format is guessed from the file extension if --output is not specified",
    )
    parser.add_argument(
        "-c",
        "--config",
        help="configuration file for generator; configuration is Python code which "
        "interacts with the variable `model` that represents the model instance",
    )

    args = parser.parse_args()

    if args.version:
        print(f"chromo {version}")
        raise SystemExit

    if args.seed <= 0:
        args.seed = int.from_bytes(os.urandom(4), "little")

    if args.number <= 0:
        raise SystemExit("Error: number must be positive")

    try:
        model_number = int(args.model)
        Model = MODELS[model_number]
    except KeyError:
        raise SystemExit(f"Error: model={args.model} is invalid ({VALID_MODELS})")
    except ValueError:
        # args.model is not a number.
        # Find model that matches string spec.
        matches = []
        x = args.model.lower()
        for M in get_all_models():
            y = M.label.lower()
            # If exact match found, stop search. This rule is necessary to e.g. select
            # Sibyll-2.3, since Sibyll-2.3d would also match otherwise, always leading to
            # ambiguity.
            if x == y:
                Model = M
                break
            if tolerant_string_match(x, y):
                matches.append(M)
        else:
            if len(matches) == 1:
                Model = matches[0]
            elif len(matches) == 0:
                raise SystemExit(
                    f"Error: model={args.model} has no match ({VALID_MODELS})"
                )
            else:
                raise SystemExit(
                    f"Error: model={args.model} is ambiguous, "
                    f"matches ({', '.join(v.label for v in matches)})"
                )
    args.model = Model

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
        raise SystemExit(
            "Error: you need to specify particle momenta or the CMS energy"
        )

    if pr != 0 or ta != 0:
        # compute sqrt(s)
        m1 = Particle.from_pdgid(args.projectile_id).mass * MeV
        m2 = Particle.from_pdgid(args.target_id).mass * MeV
        e1 = sqrt(m1**2 + pr**2)
        e2 = sqrt(m2**2 + ta**2)
        a = e1 + e2
        b = pr + ta
        s = (a + b) * (a - b)
        if s <= 0:
            raise SystemExit("Error: s <= 0")
        args.sqrts = sqrt(s)

    if args.out == "-":
        args.output = "hepmc"
        raise SystemExit("Error: output to stdout not yet supported")
    else:
        if args.out:  # filename was provided
            args.out = Path(args.out)
            # try to get format from filename extension
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
            mn = Model.pyname.lower()
            fn = f"chromo_{mn}_{args.seed}_{int(pid1)}_{int(pid2)}_{en}.{ext}"
            odir = os.environ.get("CRMC_OUT", ".")
            args.out = Path(odir) / fn
            # append extension
            if not args.out.suffixes[-2:]:
                ext = "." + args.output
                if ext.endswith("gz"):
                    ext = ext[:-2] + ".gz"
                args.out = Path(args.out).with_suffix(ext)

    if args.output not in FORMATS:
        raise SystemExit(f"Error: unknown format {args.output} ({VALID_FORMATS})")

    configuration = ""
    if args.config:
        fn = Path(args.config)
        if not fn.exists():
            raise SystemExit(f"Error: configuration file {args.config} does not exist")

        lines = ["def configure(model):\n"]
        for line in open(fn):
            lines.append(f"    {line}")
        configuration = "".join(lines)

    return args, configuration


def main():
    from rich.console import Console
    from rich.panel import Panel

    args, configuration = parse_arguments()

    p1 = Particle.from_pdgid(args.projectile_id)
    p2 = Particle.from_pdgid(args.target_id)
    pr = args.projectile_momentum
    ta = args.target_momentum

    msg = f"""\
  [repr.str]Model[/repr.str]\t\t[bold]{args.model.label}[/bold]
  [repr.str]Projectile[/repr.str]\t[bold]{p1.name}[/bold]\
 ([repr.number]{int(args.projectile_id)}[/repr.number])
  [repr.str]Target[/repr.str]\t[bold]{p2.name}[/bold]\
 ([repr.number]{int(args.target_id)}[/repr.number])\
"""
    if pr != 0 or ta != 0:
        msg += f"""
  [repr.str]Projectile momentum[/repr.str]\t[magenta]{pr:g} GeV/c[/magenta]
  [repr.str]Target momentum[/repr.str]\t[magenta]{ta:g} GeV/c[/magenta]\
"""
    msg += f"""
  [repr.str]sqrt(s)[/repr.str]\t[magenta]{args.sqrts:g} GeV[/magenta]
  [repr.str]Collisions[/repr.str]\t[repr.number]{args.number}[/repr.number]
  [repr.str]Seed[/repr.str]\t\t[repr.number]{args.seed}[/repr.number]
  [repr.str]Format[/repr.str]\t{args.output}\
"""
    if args.config:
        msg += f"\n  [repr.str]Configuration[/repr.str]\t{args.config}"

    console = Console()
    console.print(
        Panel(msg, title=f"[bold]chromo [green]{version}[/green][/bold]", width=78)
    )

    if pr > 0 and ta == 0:  # fixed target mode
        evt_kin = FixedTarget(Momentum(pr), args.projectile_id, args.target_id)
    else:  # cms mode
        evt_kin = CenterOfMass(args.sqrts, args.projectile_id, args.target_id)

    task_id = None
    try:
        model = args.model(evt_kin, seed=args.seed)
        if configuration:
            try:
                d = {}
                exec(configuration, d)
                d["configure"](model)
            except Exception:
                print(f"Error in configuration code:\n\n{configuration}")
                raise
        ofile = FORMATS[args.output](args.out, model)
        with ofile:
            # workaround: several models generate extra print when first
            # event is generated, this interferes with progress bar so we
            # create bar only after second event is generated
            with Progress(
                MofNCompleteColumn(),
                BarColumn(),
                TaskProgressColumn(),
                "ETA",
                TimeRemainingColumn(elapsed_when_finished=True),
                SpeedColumn(),
            ) as bar:
                for event in model(args.number):
                    ofile.write(event)
                    if task_id is None:
                        task_id = bar.add_task("", total=args.number)
                    bar.advance(task_id, 1)
    except Exception:
        if int(os.environ.get("DEBUG", "0")) > 0:
            raise

        import traceback

        msg = traceback.format_exc(chain=False, limit=0)
        raise SystemExit(msg.strip())


if __name__ == "__main__":
    main()
