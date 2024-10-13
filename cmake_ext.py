"""
Build an extension by calling cmake.

The actual build logic is in the CMakeLists.txt.

The following environment variables change behavior.

DEBUG=1                             : Compile in debug mode.
VERBOSE=1                           : Let cmake print commands it is calling,
                                      useful to debug the compiler options.
CMAKE_GENERATOR=<Generator>         : Specify which generator to use.
CMAKE_ARGS=<cmake args>             : Pass additional cmake arguments.
CMAKE_BUILD_PARALLEL_LEVEL=<number> : Compile in parallel with number threads.
                                      Default is to use number of CPU cores.
CMAKE_EXECUTABLE=<path>             : Path to cmake executable
"""

import os
import re
import subprocess as subp
import sys
from pathlib import Path
import sysconfig as sc
import platform

from setuptools import Extension
from setuptools.command.build_ext import build_ext

cwd = Path(__file__).parent


# Normal print does not work while CMakeBuild is running
def force_print(msg):
    subp.run(["echo", msg])


def cache_value(key, s):
    m = re.search(key + r":[A-Z]+=(.*)" + "$", s, flags=re.MULTILINE)
    assert m, f"{key} is not a cached cmake variable"
    return m.group(1)


def get_models():
    # For convenience, support building other models via models.cfg.
    # models.cfg is not tracked by git, so can be freely modified.
    # If models.cfg exists, it overrides the standard targets to
    # compile. It can be used to compile non-standard models or
    # compile only a subset of all models.
    #
    # models.cfg example:
    # -----
    # sib23c00
    # sib23c02
    # sib23c03
    # dev_dpmjetIII193=/full/path/to/dir/dpmjetIII-19.3
    # dev_sib23d=/full/path/to/dir/sibyll
    # ----
    models = {}

    for fn in (cwd / "models.cfg", cwd / "default_models.cfg"):
        if not fn.exists():
            continue
        with open(fn) as f:
            for model in f:
                model = model.strip()
                if not model or model.startswith("#"):
                    continue
                model, *mpath = model.split("=")
                models[model] = mpath[0] if mpath else None
        break

    # urqmd34, epos and pythia8 don't build correctly on Windows
    if platform.system() == "Windows":
        for model in list(models.keys()):
            if any(
                [
                    model.lower().strip().startswith(mstr)
                    for mstr in [
                        "urqmd34",
                        "eposlhc",
                        "pythia8",
                        "dpmjetIII19",
                        "phojet19",
                    ]
                ]
            ):
                del models[model]

    return models


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = Path(sourcedir).absolute()


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = Path(self.get_ext_fullpath(ext.name)).parent.absolute()

        cmake_exe = os.environ.get("CMAKE_EXECUTABLE", "cmake")

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        verbose = int(os.environ.get("VERBOSE", 0))

        cfg = "Debug" if debug else "Release"

        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        cmake_args = []
        if cmake_generator:
            cmake_args.append("-G" + cmake_generator)

        cmake_args += [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}/",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
        ]

        models = get_models()
        for model, path in models.items():
            if path is not None:
                cmake_args.append(f"-DBUILD_{model}={Path(path)}")

        # Arbitrary CMake arguments added via environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        build_args = ["--config", cfg]  # needed by some generators, e.g. on Windows

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += [f"-DCMAKE_OSX_ARCHITECTURES={';'.join(archs)}"]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators. CMake 3.12+ only. If not set, use number of
        # CPU cores.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            njobs = os.cpu_count() or 1
            build_args += [f"-j{njobs}"]

        build_temp = Path(self.build_temp)
        if not build_temp.exists():
            build_temp.mkdir(parents=True)

        # cmake setup must be run again if external settings change
        cmake_cache = build_temp / "CMakeCache.txt"
        if cmake_cache.exists():
            with cmake_cache.open() as f:
                s = f.read()
            cached_generator = cache_value("CMAKE_GENERATOR", s)
            cached_cfg = cache_value("CMAKE_BUILD_TYPE", s)
            disagreement = [
                (cmake_generator and cached_generator != cmake_generator),
                cached_cfg != cfg,
            ]
            for arg in cmake_args:
                if arg.startswith("-D"):
                    key, value = arg[2:].split("=")
                    cached_value = cache_value(key, s)
                    # Change \ to / in case of Windows
                    disagreement.append(
                        Path(value).absolute() != Path(cached_value).absolute()
                    )
            if any(disagreement):
                cmake_cache.unlink()

        # run cmake setup only once
        if not cmake_cache.exists():
            cmd = [cmake_exe, str(ext.sourcedir)] + cmake_args
            if verbose:
                force_print(" ".join(cmd))
            r = subp.run(cmd, cwd=build_temp)
            if r.returncode != 0 and cmake_cache.exists():
                cmake_cache.unlink()
                raise SystemExit("Error: CMake configuration failed")

        # This is a hack to make the parallel build faster.
        # Compile all targets on first call to make better use of parallization.
        # Skip output that already exists, but run at least once since only
        # cmake can detect which existing targets need to be recompiled.
        target = ext.name.split(".")[-1]
        suffix = sc.get_config_var("EXT_SUFFIX") or sc.get_config_var("SO")
        output_file = Path(extdir) / (target + suffix)
        targets = [f"_{m}" for m in models]
        if not output_file.exists() or target == targets[0]:  # any target is fine
            cmd = [cmake_exe, "--build", ".", "-t"] + targets + build_args
            if verbose:
                force_print(" ".join(cmd))
            r = subp.run(cmd, cwd=build_temp)
            if r.returncode != 0 and cmake_cache.exists():
                cmake_cache.unlink()
                raise SystemExit("Error: CMake configuration is faulty")
