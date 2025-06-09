#!/usr/bin/env python3
"""Script to list EPOS-LHC source files for the chromo build system."""

import argparse
import glob
import re
from pathlib import Path


def list_files(directory):
    """List EPOS-LHC source files in a directory, filtering out excluded patterns."""
    fdir = Path(directory)
    if not fdir.is_dir():
        msg = f"Directory {fdir} does not exist or is not a directory."
        raise ValueError(msg)

    sources = []
    # Collect all .f files from the sources subdirectory
    sources.extend(glob.glob(str(fdir / "sources/*.f")))

    # Filter out excluded files
    exclude_patterns = [
        r"epos_example\.f",  # Example file
        r"epos-random\.f",  # Random number generator (replaced by dummy)
    ]

    filtered_sources = []
    for src in sources:
        excluded = False
        for pattern in exclude_patterns:
            if re.search(pattern, src):
                excluded = True
                break
        if not excluded:
            filtered_sources.append(src)

    # Add the dummy random number generator
    filtered_sources.append(str(fdir / "epos-random-dummy.f"))

    # Sort for consistent output
    filtered_sources.sort()

    for src in filtered_sources:
        print(src)  # noqa: T201


def list_interface_files(directory):
    """List EPOS-LHC interface source files specifically used for f2py interface."""
    fdir = Path(directory)
    if not fdir.is_dir():
        msg = f"Directory {fdir} does not exist or is not a directory."
        raise ValueError(msg)

    # Specific files used for f2py interface
    interface_files = [
        fdir / "sources/epos-bas-lhc.f",
        fdir / "sources/epos-ids-lhc.f",
        fdir / "sources/epos_interface.f",
        fdir / "epos-random-dummy.f",
    ]

    # Only include files that actually exist
    existing_files = [str(f) for f in interface_files if f.exists()]
    existing_files.sort()

    for src in existing_files:
        print(src)  # noqa: T201


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="List EPOS-LHC source files.")
    parser.add_argument(
        "--directory",
        type=str,
        required=True,
        help="Directory containing the EPOS-LHC source files.",
    )
    parser.add_argument(
        "--interface-only",
        action="store_true",
        help="List only interface files used for f2py compilation.",
    )

    args = parser.parse_args()

    if args.interface_only:
        list_interface_files(args.directory)
    else:
        list_files(args.directory)
