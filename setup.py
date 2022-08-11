from pathlib import Path
from setuptools import setup
import sys

cwd = Path(__file__).parent

sys.path.append(str(cwd))
from cmake_ext import CMakeExtension, CMakeBuild  # noqa: E402

setup(
    zip_safe=False,
    ext_modules=[
        CMakeExtension("impy.models.eposlhc"),
        CMakeExtension("impy.models.sib21"),
        CMakeExtension("impy.models.sib23d"),
    ],
    cmdclass={"build_ext": CMakeBuild},
)
