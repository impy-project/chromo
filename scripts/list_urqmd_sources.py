#!/usr/bin/env python3
"""
Generate URQMD source-file lists for Meson.

 * With no flag  → prints **all** .f/.f90/.fpp sources that belong
   in the static library (minus the *global* exclude list).

 * With --interface-only → prints only those files that should be
   handed to f2py when building the Python front-end.

The file names and rules are taken 1-to-1 from the explicit lists
that used to live in meson.build (urqmd_sources_f / urqmd_sources_f90 /
urqmd_ignore_in_interface_sources).
"""

import argparse
from pathlib import Path

# ---  hard-coded lists copied from meson.build  ------------------------------------
SOURCES_F = [
    "1fluid.f",
    "bessel.f",
    "delpart.f",
    "getmass.f",
    "hepcmp.f",
    "iso.f",
    "numrec.f",
    "pythia6409.f",
    "siglookup.f",
    "upmerge.f",
    "addpart.f",
    "blockres.f",
    "coload.f",
    "detbal.f",
    "getspin.f",
    "hepnam.f",
    "ityp2pdg.f",
    "output.f",
    "string.f",
    "urqmd.f",
    "angdis.f",
    "boxprg.f",
    "dwidth.f",
    "init.f",
    "jdecay2.f",
    "paulibl.f",
    "saveinfo.f",
    "tabinit.f",
    "whichres.f",
    "anndec.f",
    "cascinit.f",
    "dectim.f",
    "error.f",
    "hepchg.f",
    "input.f",
    "make22.f",
    "proppot.f",
    "scatter.f",
    "uhmerge.f",
    "urqinit.f",
]
SOURCES_F90 = ["CFmax.f90", "quadri.f90", "cornelius.f90"]
CHROMO_SRC = ["chromo_urqmd.f"]  # always included at the end

# files to exclude from the *interface* list (they still compile into the lib)
INTERFACE_EXCLUDE = {"newpart.f", "uhmerge.f", "iso.f"}

# -----------------------------------------------------------------------------------


def collect(root: Path, interface_only: bool):
    # Build the list in the exact order used previously
    seq = [
        *(root / f for f in SOURCES_F),
        *(root / f for f in SOURCES_F90),
        *(root.parent / f for f in CHROMO_SRC),
    ]  # chromo_urqmd sits one level up

    for p in seq:
        if interface_only and p.name in INTERFACE_EXCLUDE:
            continue
        print(p.as_posix())  # noqa: T201


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--directory", required=True, help="Path to urqmd-3.4/sources")
    ap.add_argument(
        "--interface-only",
        action="store_true",
        help="Emit only the subset handed to f2py",
    )
    args = ap.parse_args()

    collect(Path(args.directory).resolve(), args.interface_only)


if __name__ == "__main__":
    main()
