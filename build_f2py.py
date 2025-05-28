#!/usr/bin/env python
"""Build a Fortran extension using f2py with the meson backend."""
from __future__ import annotations
import sys
import subprocess
import shutil
from pathlib import Path


def main() -> None:
    if len(sys.argv) < 3:
        print("Usage: build_f2py.py <output> <sources...>")
        sys.exit(1)
    output = Path(sys.argv[1]).resolve()
    modulename = output.stem.split(".")[0]
    sources = [str(Path(s)) for s in sys.argv[2:]]

    fcflags = (
        "-std=legacy -fno-second-underscore -Wno-uninitialized "
        "-Wno-argument-mismatch -fallow-argument-mismatch -cpp"
    )
    if shutil.which("mold"):
        fcflags += " -fuse-ld=mold"

    command = [
        sys.executable,
        "-m",
        "numpy.f2py",
        "--backend",
        "meson",
        f"--f90flags={fcflags}",
        f"--f77flags={fcflags}",
        "-c",
        "-m",
        modulename,
    ] + sources

    log_file = output.with_suffix(".log")
    with log_file.open("w") as log:
        proc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise SystemExit(proc.returncode)

    built = None
    for pattern in (f"{modulename}*.so", f"{modulename}*.pyd"):
        matches = list(Path(".").glob(pattern))
        if matches:
            built = matches[0]
            break
    if not built:
        raise FileNotFoundError(f"Built module {modulename} not found")
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(built), output)
    print(f"Built {output}")


if __name__ == "__main__":
    main()
