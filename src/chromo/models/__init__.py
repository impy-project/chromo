import platform
from chromo.models.sophia import Sophia20
from chromo.models.sibyll import (
    Sibyll21,
    Sibyll23,
    Sibyll23c,
    Sibyll23d,
)

from chromo.models.dpmjetIII import DpmjetIII306
from chromo.models.epos import EposLHC
from chromo.models.qgsjet import QGSJet01d, QGSJetII03, QGSJetII04
from chromo.models.dpmjetIII import DpmjetIII191, DpmjetIII193
from chromo.models.phojet import Phojet112, Phojet191, Phojet193
from chromo.models.urqmd import UrQMD34
from chromo.models.pythia6 import Pythia6
from chromo.models.pythia8 import Pythia8

__all__ = (
    "Sophia20",
    "Sibyll21",
    "Sibyll23",
    "Sibyll23c",
    "Sibyll23d",
    "DpmjetIII306",
    "EposLHC",
    "QGSJet01d",
    "QGSJetII03",
    "QGSJetII04",
    "DpmjetIII191",
    "DpmjetIII193",
    "Phojet112",
    "Phojet191",
    "Phojet193",
    "UrQMD34",
    "Pythia6",
    "Pythia8",
)


# This is a workaround of a problem with "fixed" wheels
# which is produced by "delocate" package on cibuildwheels
# It changes Python -> the path of libpython in the "*.so" files
# using standard "otool" and "install_name_tool" routines
def fix_macos_installation(logfile):
    import subprocess
    import re
    from pathlib import Path
    from distutils.sysconfig import get_config_var

    models_dir = Path(__file__).parent
    for ext_file in models_dir.glob("*.so"):
        ext_file = str(ext_file)

        # Get libraries required by extension module
        cmd = ["otool", "-L", ext_file]
        result = subprocess.run(cmd, stdout=subprocess.PIPE)
        lib_deps = result.stdout.decode("utf-8")

        message = f"\n----------\nFixing {ext_file}:\nBefore\n{lib_deps}"
        with open(logfile, "a") as lf:
            lf.write(message)

        # Search for python....dylib or Python shared library
        search_res = re.search(
            "^(.*python.*dylib|.*Python).*$", lib_deps, flags=re.MULTILINE
        )

        if search_res is None:
            message = (
                'NOT FIXED: Paths  with "Python" or "python*.dylib" strings '
                "isn't found inside shared library file"
            )
            with open(logfile, "a") as lf:
                lf.write(message)
            continue

        python_lib_in_file = search_res.group(1).strip()
        python_lib_real = str(next(Path(get_config_var("LIBDIR")).glob("libpython*")))

        # Change paths
        cmd = [
            "install_name_tool",
            "-change",
            python_lib_in_file,
            python_lib_real,
            ext_file,
        ]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Since Big Sur, codesigning for arm64 is much more strict than it is for x64.
        # https://stackoverflow.com/a/71753248
        cmd = [
            "codesign",
            "--force",
            "-s",
            "-",
            ext_file,
        ]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        cmd = ["otool", "-L", ext_file]
        result = subprocess.run(cmd, stdout=subprocess.PIPE)
        lib_deps = result.stdout.decode("utf-8")

        message = f"\nAfter\n{lib_deps}"
        with open(logfile, "a") as lf:
            lf.write(message)


if platform.system() == "Darwin":
    from pathlib import Path

    logfile = Path(__file__).parent / "macos_fixed_ext.log"
    if not logfile.exists():
        fix_macos_installation(logfile)
