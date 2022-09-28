# The module is for development purposes
# It is used for monitoring environment variables
# which turn on additional settings:
#
# IMPY_DEV_EXTRA_MODELS=1 pip install --prefer-binary\
# --no-build-isolation --no-deps -v -e .
# builds additional models
#
# IMPY_DEV_DPMJETIII193=/full/path/to/dpmjetIII193 \
# pip install --prefer-binary --no-build-isolation --no-deps -v -e .
# builds dpmjetIII193 from given directory
#
# IMPY_DEV_PYF_GENERATION pip install --prefer-binary\
# --no-build-isolation --no-deps -v -e .
# generates *.pyf files during building
#
from pathlib import Path
import sys
import os
import json

cwd = Path(__file__).parent


def development_settings(models):
    cache_file = "env_variables_cache.json"

    # Monitor specific environment variables and update
    # file "env_variables_cache.json" if they change
    # Define the environment variables which we monitor

    env_variables = {
        "IMPY_DEV_EXTRA_MODELS": None,
        "IMPY_DEV_DPMJETIII193": None,
    }

    # Read them if already exists
    if Path(cache_file).exists():
        with open(cache_file, "r") as file:
            for item in json.load(file).items():
                env_variables[item[0]] = item[1]
    else:
        with open(cache_file, "w") as file:
            json.dump(env_variables, file)

    # Check if the variables are updated
    env_updated = False
    for evar in env_variables.items():
        if evar[0] in os.environ:
            if evar[1] != os.environ[evar[0]]:
                env_updated = True
                env_variables[evar[0]] = os.environ[evar[0]]
        else:
            if evar[1]:
                env_updated = True
                env_variables[evar[0]] = None

    # Write to the file if the variables are updated
    if env_updated:
        with open(cache_file, "w") as file:
            json.dump(env_variables, file)

    if env_variables["IMPY_DEV_DPMJETIII193"]:
        models.append("dev_dpmjetIII193")

    if env_variables["IMPY_DEV_EXTRA_MODELS"]:
        models.append("sib23c00")
        models.append("sib23c02")
        models.append("sib23c03")

    # Check for existence of f2py files
    pyf_files = []
    for model in models:
        pyf_files.append(f"_{model}.pyf")
        pyf_files.append(f"_{model}module.c")
        pyf_files.append(f"_{model}-f2pywrappers.f")

    f2py_dir = f"{str(cwd)}/src/f2py"
    # Touch file if not all f2py files in the f2py directory "f2py_dir"
    # lfile = f"{str(cwd)}/generate_f2py_files.dat"
    # info_string = "This file is for automatic regeneration of f2py files by CMake"

    # if not Path(lfile).exists():
    #     with open(lfile, "w") as f:
    #         f.write(info_string)
    # else:
    flist = [str(f.name) for f in Path(f2py_dir).iterdir()]
    for filename in pyf_files:
        if filename not in flist:
            print(f"{filename} is not in flist")
            print("TOUCHHHHHHHHHHH")
            Path(cache_file).touch()
            break
            # with open(lfile, "w") as f:
            #     f.write(info_string)
            # break
