from pathlib import Path
from setuptools import setup
import sys
import subprocess as subp
import json

cwd = Path(__file__).parent

sys.path.append(str(cwd))
from cmake_ext import CMakeExtension, CMakeBuild  # noqa: E402

# make sure that submodules are up-to-date,
# it is a common error to forget this when
# switching between development branches
subp.check_call(["git", "submodule", "update"])

models = [
    "eposlhc",
    "sib21",
    "sib23",
    "sib23d",
    "sib23c01",
    "qgs01",
    "qgsII03",
    "qgsII04",
    "urqmd34",
    "pythia6",
    "sophia",
    "dpmjet306",
    "phojet112",
    "phojet191",
    "phojet193",
    "dpmjetIII191",
    "dpmjetIII193",
]

# for convenience, support building extra models via extra.cfg
# extra.cfg is not tracked by git, so can be freely modified
extra_cfg = cwd / "extra.cfg"

if extra_cfg.exists():
    with open(extra_cfg) as f:
        for line in f:
            models.append(line.strip())

ext_modules = []
for model in models:
    ext_modules.append(CMakeExtension(f"impy.models._{model}"))

setup(
    zip_safe=False,
    ext_modules=ext_modules,
    cmdclass={"build_ext": CMakeBuild},
)
