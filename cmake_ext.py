import os
import re
import subprocess
import sys
from pathlib import Path

from setuptools import Extension
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = Path(sourcedir).absolute()


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = Path(self.get_ext_fullpath(ext.name)).parent.absolute()

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}/",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
        ]

        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        build_args = ["--config", cfg]  # needed by some generators, e.g. on Windows

        if self.compiler.compiler_type == "msvc":
            # CMake allows an arch-in-generator style for backward compatibility
            contains_arch = any(x in cmake_generator for x in ("ARM", "Win64"))

            # Specify the arch if using MSVC generator, but only if it doesn't
            # contain a backward-compatibility arch spec already in the
            # generator name.
            if not contains_arch:
                # Convert distutils Windows platform specifiers to CMake -A arguments
                arch = {
                    "win32": "Win32",
                    "win-amd64": "x64",
                    "win-arm32": "ARM",
                    "win-arm64": "ARM64",
                }[self.plat_name]
                cmake_args += ["-A", arch]

            cmake_args += [f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"]

        elif sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += [f"-DCMAKE_OSX_ARCHITECTURES={';'.join(archs)}"]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # CMake 3.12+ only.
            njobs = self.parallel or os.cpu_count() or 1
            build_args += [f"-j{njobs}"]

        build_temp = Path(self.build_temp)
        if not build_temp.exists():
            build_temp.mkdir(parents=True)

        # cmake setup must run again if python path has changed
        cmake_cache = build_temp / "CMakeCache.txt"
        if cmake_cache.exists():
            with cmake_cache.open() as f:
                m = re.search("PYTHON_EXECUTABLE:FILEPATH=([^\s]+)", f.read())
                cached_python_path = m.group(1)
                if cached_python_path != sys.executable:
                    cmake_cache.unlink()
        if not cmake_cache.exists():
            print(f"cmake args: {' '.join(cmake_args)}")
            print(f"build args: {' '.join(build_args)}")

            subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=build_temp)

        target = ext.name.split(".")[-1]

        subprocess.check_call(
            ["cmake", "--build", ".", "--target", target] + build_args,
            cwd=build_temp,
        )
