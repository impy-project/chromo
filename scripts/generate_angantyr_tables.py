"""Generate Angantyr precomputed initialization tables.

Replicates the logic of Pythia8 examples/main424.cc (option 2) to produce
an InitDefaultAngantyr.cmnd file for arbitrary CMS energy ranges.  The
bundled file covers 20 GeV – 20 PeV; use this script if you need tables
at higher energies (e.g. beyond 200 EeV lab / ~600 TeV CMS).

This initialization run is slow (hours for large grids).  The output file
can be placed in the setups/ directory to replace the bundled tables.

Usage:
    python scripts/generate_angantyr_tables.py [options]

Options:
    --ecm-max FLOAT   Maximum CMS energy in GeV (default: 1e8 = 100 PeV)
    --ecm-min FLOAT   Minimum CMS energy in GeV (default: 10.0)
    --grid-pts INT    Number of grid points for sigma fitting (default: 25)
    --output FILE     Output .cmnd file path (default: InitDefaultAngantyr_custom.cmnd)
    --mpi-init FILE   Path to MPI init file (default: bundled InitDefaultMPI.cmnd)
"""

import argparse
import sys
from os import environ
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


def find_setup_file(filename):
    """Locate a Pythia8 setup file in the source tree."""
    from chromo.models.pythia8 import Pythia8
    from chromo.util import _cached_data_dir

    datdir = Path(_cached_data_dir(Pythia8._data_url) + "xmldoc")
    candidates = [
        datdir.parent / "setups" / filename,
        Path(__file__).parent.parent
        / "src/cpp/pythia83/share/Pythia8/setups"
        / filename,
    ]
    for p in candidates:
        if p.exists():
            return str(p)
    raise FileNotFoundError(
        f"Could not find {filename}. "
        "Ensure the Pythia8 submodule is checked out or add setups/ to the data bundle."
    )


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--ecm-max",
        type=float,
        default=1e8,
        help="Maximum CMS energy in GeV (default: 1e8)",
    )
    parser.add_argument(
        "--ecm-min",
        type=float,
        default=10.0,
        help="Minimum CMS energy in GeV (default: 10.0)",
    )
    parser.add_argument(
        "--grid-pts",
        type=int,
        default=25,
        help="Number of sigma-fit grid points (default: 25)",
    )
    parser.add_argument(
        "--output",
        default="InitDefaultAngantyr_custom.cmnd",
        help="Output file path (default: InitDefaultAngantyr_custom.cmnd)",
    )
    parser.add_argument(
        "--mpi-init", default="", help="Path to MPI init .cmnd file (default: bundled)"
    )
    args = parser.parse_args()

    from chromo.models.pythia8 import Pythia8
    from chromo.util import _cached_data_dir

    datdir = _cached_data_dir(Pythia8._data_url) + "xmldoc"

    if "PYTHIA8DATA" in environ:
        del environ["PYTHIA8DATA"]

    import _pythia8

    print(f"Initializing Pythia8 with datdir={datdir}")  # noqa: T201
    pythia = _pythia8.Pythia(datdir, True)

    mpi_init = args.mpi_init or find_setup_file("InitDefaultMPI.cmnd")
    print(f"Using MPI init file: {mpi_init}")  # noqa: T201

    config = [
        f"include = {mpi_init}",
        "Beams:allowIDAswitch = on",
        "Beams:allowVariableEnergy = on",
        f"Beams:eCM = {args.ecm_max:.4e}",
        "HeavyIon:SasdMpiReuseInit = -1",
        "HeavyIon:SigFitReuseInit = -1",
        f"HeavyIon:varECMMax = {args.ecm_max:.4e}",
        f"HeavyIon:varECMMin = {args.ecm_min:.4e}",
        f"HeavyIon:varECMSigFitNPts = {args.grid_pts}",
        "HeavyIon:varECMStepwiseEvolve = 2",
        "MultipartonInteractions:nSample = 1000000",
    ]

    for line in config:
        if not pythia.readString(line):
            print(
                f"WARNING: readString({line!r}) returned False", file=sys.stderr
            )

    print("Running Pythia8 initialization (this may take hours) ...")  # noqa: T201
    if not pythia.init():
        print("ERROR: Pythia8 initialization failed", file=sys.stderr)  # noqa: T201
        sys.exit(1)

    print(f"Writing output to {args.output!r} ...")  # noqa: T201
    if not pythia.settings.writeFile(args.output, True):
        print(
            f"ERROR: writeFile({args.output!r}) failed", file=sys.stderr
        )
        sys.exit(1)

    print(f"Done. Tables written to {args.output}")  # noqa: T201
    print(  # noqa: T201
        "Place this file in the setups/ directory and update AngantyrCascade.cmnd to include it."
    )


if __name__ == "__main__":
    main()
