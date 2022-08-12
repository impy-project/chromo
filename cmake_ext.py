"""
Build an extension by calling cmake.

The actual build logic is in the CMakeLists.txt.

The following environment variables change behavior.

DEBUG=1                             : Compile in debug mode.
VERBOSE=1                           : Let cmake print the commands it is calling,
                                      useful to debug the compiler options.
CMAKE_GENERATOR=<Generator>         : Specify which generator to use.
CMAKE_ARGS=<cmake args>             : Pass additional cmake arguments.
CMAKE_BUILD_PARALLEL_LEVEL=<number> : Compile in parallel with number threads.
"""

import os
import re
import subprocess as subp
import sys
from pathlib import Path
import sysconfig as sc

from setuptools import Extension
from setuptools.command.build_ext import build_ext


# Use Ninja if it is available instead of cmake default generator
try:
    out = subp.check_output(["ninja", "--version"])
    print(f"Ninja generator detected: {out.decode().strip()}")
    default_cmake_generator = "Ninja"
except Exception:
    default_cmake_generator = ""


# Normal print does not work while CMakeBuild is running
def force_print(msg):
    subp.run(["echo", msg])


def cache_value(key, s):
    m = re.search(key + r":[A-Z]+=([^\s]*)", s)
    assert m
    return m.group(1)


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = Path(sourcedir).absolute()


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = Path(self.get_ext_fullpath(ext.name)).parent.absolute()

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        verbose = int(os.environ.get("VERBOSE", 0))

        cfg = "Debug" if debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", default_cmake_generator)

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}/",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
        ]

        if cmake_generator:
            cmake_args.append("-G" + cmake_generator)

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

        # cmake setup must be run again if external settings have changed
        cmake_cache = build_temp / "CMakeCache.txt"
        if cmake_cache.exists():
            with cmake_cache.open() as f:
                s = f.read()
                cached_python_path = cache_value("PYTHON_EXECUTABLE", s)
                cached_generator = cache_value("CMAKE_GENERATOR", s)
                cached_cfg = cache_value("CMAKE_BUILD_TYPE", s)
                if (
                    cached_python_path != sys.executable
                    or (cmake_generator and cached_generator != cmake_generator)
                    or cached_cfg != cfg
                ):
                    cmake_cache.unlink()

        # run cmake setup only once
        if not cmake_cache.exists():
            cmd = ["cmake", str(ext.sourcedir)] + cmake_args
            if verbose:
                force_print(" ".join(cmd))
            subp.check_call(cmd, cwd=build_temp)

        # This is a hack to make the parallel build faster.
        # Compile all targets on first call to make better use of parallization.
        # Skip if output that already exists, but run at least once since only
        # cmake can detect which existing targets need to be recompiled.
        target = ext.name.split(".")[-1]
        suffix = sc.get_config_var("EXT_SUFFIX") or sc.get_config_var("SO")
        output_file = Path(extdir) / (target + suffix)
        if not output_file.exists() or target == "_eposlhc":  # any one target is fine
            cmd = ["cmake", "--build", "."] + build_args
            if verbose:
                force_print(" ".join(cmd))
            subp.check_call(cmd, cwd=build_temp)
