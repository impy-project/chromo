import os
import sys
from importlib.metadata import version

# On Windows, add the MinGW runtime DLL directory so that gfortran-compiled
# extension modules (.pyd) can find libgfortran-5.dll and libgcc_s_seh-1.dll.
if sys.platform == "win32" and hasattr(os, "add_dll_directory"):
    import shutil

    _gfortran = shutil.which("gfortran")
    if _gfortran:
        os.add_dll_directory(os.path.dirname(_gfortran))

from chromo import constants, kinematics, models

__version__ = version("chromo")

debug_level = int(os.environ.get("DEBUG", "0"))

__all__ = ["__version__", "constants", "debug_level", "kinematics", "models"]
