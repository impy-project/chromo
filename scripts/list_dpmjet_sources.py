import argparse
import glob
import re
from pathlib import Path


def list_files(directory):
    """List files in a directory matching a given pattern."""
    fdir = Path(directory)  # "src/fortran/dpmjetIII-19.3")
    if not fdir.is_dir():
        msg = f"Directory {fdir} does not exist or is not a directory."
        raise ValueError(msg)
    sources = []
    sources.extend(glob.glob(str(fdir / "src/phojet/*.f")))
    sources.extend(glob.glob(str(fdir / "src/pythia/*.f")))
    sources.extend(glob.glob(str(fdir / "src/dpmjet/*.f")))
    sources.append(str(fdir / "common/dummies.f"))

    # Filter out excluded files
    exclude_files = [
        r"DT_RNDM\.f",
        r"DT_RNDMST\.f",
        r"DT_RNDMTE\.f",
        r"PYR\.f",
    ]
    filtered_sources = []
    for src in sources:
        excluded = False
        for pattern in exclude_files:
            if re.search(pattern, src):
                excluded = True
                break
        if not excluded:
            filtered_sources.append(src)

    # Add dummies.f
    filtered_sources.append(str(fdir / "common/dummies.f"))

    for src in filtered_sources:
        print(src)  # noqa: T201


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="List DPMJET source files.")
    parser.add_argument(
        "--directory",
        type=str,
        help="Directory containing the DPMJET source files.",
    )

    args = parser.parse_args()

    list_files(args.directory)
