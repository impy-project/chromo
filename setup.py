from pathlib import Path
from setuptools import setup
import sys
import os


cwd = Path(__file__).parent

sys.path.append(str(cwd))
from cmake_ext import CMakeExtension, CMakeBuild  # noqa: E402


ext_modules = [
    CMakeExtension("impy.models._eposlhc"),
    CMakeExtension("impy.models._sib21"),
    CMakeExtension("impy.models._sib23"),
    CMakeExtension("impy.models._sib23c00"),
    CMakeExtension("impy.models._sib23c01"),
    CMakeExtension("impy.models._sib23c02"),
    CMakeExtension("impy.models._sib23c03"),
    CMakeExtension("impy.models._sib23d"),
    CMakeExtension("impy.models._qgs01"),
    CMakeExtension("impy.models._qgsII03"),
    CMakeExtension("impy.models._qgsII04"),
    CMakeExtension("impy.models._urqmd34"),
    CMakeExtension("impy.models._pythia6"),
    CMakeExtension("impy.models._sophia"),
    CMakeExtension("impy.models._dpmjet306"),
    CMakeExtension("impy.models._phojet112"),
    CMakeExtension("impy.models._dpmjetIII191"),
    CMakeExtension("impy.models._dpmjetIII193"),
    CMakeExtension("impy.models._phojet191"),
    # CMakeExtension("impy.models._phojet193"),
]


# The below block is only for development

# Alwasy create file to keep value of environment variable
# DEVELOP_DPMJETIII193_SOURCE
with open("env_variables_cache.dat", "w") as file:
    file.write("")

# Check if DEVELOP_DPMJETIII193_SOURCE is defined
develop_dir = "DEVELOP_DPMJETIII193_SOURCE"
if develop_dir in os.environ:
    cached_value = ""
    if os.path.exists("env_variables_cache.dat"):
        with open("env_variables_cache.dat", "r") as file:
            cached_value = file.read()
    # Rewrite file if DEVELOP_DPMJETIII193_SOURCE has been changed
    if cached_value != os.environ[develop_dir]:
        with open("env_variables_cache.dat", "w") as file:
            file.write(os.environ[develop_dir])
    ext_modules.append(CMakeExtension("impy.models._dev_dpmjetIII193"))
# End of block for development

setup(
    zip_safe=False,
    ext_modules=ext_modules,
    cmdclass={"build_ext": CMakeBuild},
)
