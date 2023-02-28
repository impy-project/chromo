# ![](doc/chromo.svg)<br> Cosmic ray and HadROnic interactiOn MOnte-carlo frontend

This package implements are generic user interface to popular event generators used in cosmic ray and high-energy particle physics. The purpose of the package is to simplify working with simulations of particle interactions without the need to use Fortran style interfaces to event generators, 'ASCII input cards' and files or C++ dependencies.  

Simulate interactions with one of the supported event generators 

```python
import numpy as np
import chromo

# Define the parameters of the collisions
kinematics = chromo.kinematics.CenterOfMass(
    13 * chromo.constants.TeV,
    "proton", "proton")
# Create an instance of an event generator
generator = chromo.models.Sibyll23d(kinematics)

nevents = 0
average_pt = 0

# Generate 10000 events
for event in generator(10000):
    # Filter event
    event = event.final_state_charged()
    # do something with event.pid, event.eta, event.en, event.pt, etc.
    # these variables are numpy arrays, that can be histogrammed or counted like
    pt = event.pt[np.abs(event.pid) == 211]
    # The list could be empty
    if len(pt) > 0:
        nevents += 1
        average_pt += np.mean(pt)

average_pt = average_pt / nevents
print("Average pT for charged pions {0:4.3f}".format(average_pt))
```

## Supported models

- DPMJET-III 3.0.6 & PHOJET 1.12-35
- DPMJET-III 19.1 & PHOJET 19.1
- DPMJET-III 19.3 & PHOJET 19.3
- EPOS-LHC
- PYTHIA 6.4
- PYTHIA 8.3
- QGSJet-01
- QGSJet-II-03
- QGSJet-II-04
- SIBYLL-2.1
- SIBYLL-2.3
- SIBYLL-2.3c
- SIBYLL-2.3d
- SOPHIA 2.0
- UrQMD 3.4

## Installation

## Supported platforms

- Python 3.8+
- Linux, Mac OS X, Windows

### From PyPI (not yet available)

    pip install chromo

The package will be available as a pre-compiled binary wheels in the future, but for now you have to compile it from source, see next subsection.

### From source

Installation from source requires a Python installation setup for development, as well as C and Fortran compilers.

To build from source (the **recursive** flag is important to check out submodules):

    git clone --recursive https://github.com/impy-project/chromo
    cd chromo
    pip install --prefer-binary -v -e .

This takes a while. The command installs the package in editable mode (for developing the Python layer) and with verbose output, so that you can watch the compilation happening. Warnings can be ignored, but watch out for errors. Pip automatically creates a virtual environment and downloads dependencies for the build, prefering to install binary wheels of older versions of dependencies if the most recent version has no binary wheel, so you don't have to compile dependencies as well.

To run the tests or try the examples, use this modified `pip install` instead:

    pip install --prefer-binary -v -e .'[test,examples]'

which installes chromo and additional optional Python packages that are used in the tests and examples, but not required to run chromo.

#### For developers

If you want to work on the Fortran sources, it is more convenient to install a development version of Chromo with setuptools.

    python setup.py develop

You need to install the build environment manually for this to succeed. Check the key `[build-system.requires]` in `pyproject.toml` which packages are required.

Unlike the pip command, this command reuses build artefacts that were previously generated, so you don't have to recompile everything every time. Another convenience for developers is the optional file `models.cfg`. If it exists, only models are build which are listed there. See `default_models.cfg` for the full list.

#### Known issues

- On OSX
    - You need to install gcc and gfortran with homebrew.
    - Apple introduced a bug in the Xcode Command Line Tools Version 14 which produces a linker error when compiling C++ code with gcc. Until this is fixed, the workaround is to downgrade to 13.4, use this link https://download.developer.apple.com/Developer_Tools/Command_Line_Tools_for_Xcode_13.4/Command_Line_Tools_for_Xcode_13.4.dmg and turn off automatic updates in the System Settings, because otherwise your Mac will upgrade to 14 again.
- setuptools > 60 does not seem to work. Downgrade with `pip install setuptools<60` if you experience problems.

If you cannot fix the installation with these hints, please look into the subsection below which explains how to install in chromo in a verified docker environment. The docker environment has a properly set up environment verified by us, so that the installation is guaranteed to succeed.

### From source in Docker

This guide works on Linux and OSX. You need a running Docker server. Please google how to set up Docker on your machine.

    # download chromo
    git clone --recursive https://github.com/impy-project/chromo
    cd chromo

    # download linux image for x86_64 or see below
    docker pull quay.io/pypa/manylinux2014_x86_64
 
    # For aarch64 or VM on Apple Silicon use the following image and
    # replace the end of the next command accordingly.
    # docker pull quay.io/pypa/manylinux2014_aarch64
    
    # create docker instance and bind chromo directory
    docker run --rm -d -it --name chromo -v "$(pwd)":/app quay.io/pypa/manylinux2014_x86_64

    # enter your docker instance
    docker exec -it chromo /bin/bash

    cd /app

    # select python version, e.g. 3.9, and enter virtual environment
    python3.9 -m venv venv
    source venv/bin/activate

    # install chromo and dependencies (prefer binary wheels for deps)
    pip install --prefer-binary -v -e .

You can now use chromo inside the docker instance. If you run Linux, you can also make a wheel inside
docker and install it in your host.

    # inside docker
    pip install wheel
    python setup.py bdist_wheel

    # exit docker with ctrl+D
    pip install dist/*.whl

This should allow you to use chromo also outside docker. This works only if you use the same Python version inside and outside of docker.

## User interface

There are two ways to interact with the code.

1. As in the example above, via plain python in scripts or Jupyter notebooks. Look at this [example](examples/compare_models.ipynb).

2. Via a HEPMC output that can be piped in Rivet or other tools supporting the format.

## Running tests

Some notes regarding tests.

- Tests are run in parallel by default with `pytest-xdist`. To disable this, use the option `-n 0`.
- The test `test_generators` takes a long time. You can skip it with the option `-k "not test_generators"`.
- Tests which run a model do so in a separate process, because most models can only be instantiated once. This prevents using `--pdb` to start the debugger at the point of failure. You can prefix the pytest call like this `DEBUG=10 python -m pytest ...` to run the model in the current process. This will only work once for each model and lead to failures afterwards.

## Authors

- Anatoli Fedynitch
- Hans Dembinski
- Anton Prosekin
- Sonia El Hadri
- Keito Watanabe

## LICENSE

The source code of chromo is licensed under the [BSD 3-clause license (see LICENSE for detail)](LICENSE). The source codes of the event generators are individually licensed under different conditions (see the COPYING files located in the subdirectories).
