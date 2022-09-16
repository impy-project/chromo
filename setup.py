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
    CMakeExtension("impy.models._sib23d"),
    CMakeExtension("impy.models._qgs01"),
    CMakeExtension("impy.models._qgsII03"),
    CMakeExtension("impy.models._qgsII04"),
    CMakeExtension("impy.models._urqmd34"),
    CMakeExtension("impy.models._pythia6"),
    CMakeExtension("impy.models._sophia"),
    CMakeExtension("impy.models._dpmjet306"),
    CMakeExtension("impy.models._phojet112"),
    CMakeExtension("impy.models._phojet191"),
    CMakeExtension("impy.models._dpmjetIII191"),
    CMakeExtension("impy.models._dpmjetIII193"),
    # CMakeExtension("impy.models._phojet193"),
]


# The below block is only for development

# Monitor specific environment variables and update
# file "env_variables_cache.dat" if they change
# Define the environment variables which we monitor
cache_file = "env_variables_cache.dat"
env_variables = dict()
env_variables["DEVELOP_DPMJETIII193_SOURCE"] = "not defined"
env_variables["IMPY_GENERATE_PYF"] = "not defined"
env_variables["IMPY_EXTRA_MODELS"] = "not defined"

# Read them if already exists
if os.path.exists(cache_file):
    with open(cache_file, "r") as file:
        Lines = file.readlines()
        for line in Lines:
            linevals = line.split("=")
            env_variables.update({linevals[0]: linevals[1]})
else:
    with open(cache_file, "w") as file:
        for env_variable in env_variables.items():
            line_string = env_variable[0] + "=" + env_variable[1] + "\n"
            file.write(line_string)

# Check if the variables are updated
env_updated = False
for env_variable in env_variables.items():
    if env_variable[0] in os.environ:
        if env_variable[1] != os.environ[env_variable[0]]:
            env_updated = True
            env_variables.update({env_variable[0]: os.environ[env_variable[0]]})
    else:
        if env_variable[1] != "not defined":
            env_updated = True
            env_variables.update({env_variable[0]: "not defined"})

# Write to the file if the variables are updated
if env_updated:
    with open(cache_file, "w") as file:
        for env_variable in env_variables.items():
            line_string = env_variable[0] + "=" + env_variable[1] + "\n"
            file.write(line_string)

if env_variables["DEVELOP_DPMJETIII193_SOURCE"] != "not defined":
    ext_modules.append(CMakeExtension("impy.models._dev_dpmjetIII193"))


if env_variables["IMPY_EXTRA_MODELS"] != "not defined":
    ext_modules.append(CMakeExtension("impy.models._sib23c00"))
    ext_modules.append(CMakeExtension("impy.models._sib23c01"))
    ext_modules.append(CMakeExtension("impy.models._sib23c02"))
    ext_modules.append(CMakeExtension("impy.models._sib23c03"))
# End of block for development

setup(
    zip_safe=False,
    ext_modules=ext_modules,
    cmdclass={"build_ext": CMakeBuild},
)
