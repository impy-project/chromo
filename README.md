# impy - (hadronic) interaction models in python

This package implements are generic user interface to popular event generators used in cosmic ray and high-energy particle physics. The purpose of the package is to simplify working with simulations of particle interactions without the need to use Fortran style interfaces to event generators, 'ASCII input cards' and files or C++ dependencies.  

Simulate interactions with one of the supported event generators 

```python
import numpy as np
import impy

# Define the parameters of the collisions
kinematics = impy.kinematics.CenterOfMass(
    13 * impy.constants.TeV,
    "proton", "proton")
# Create an instance of an event generator
generator = impy.models.Sibyll23d(kinematics)

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

## Installation

## Supported platforms

- Python 3.6+
- Linux, Mac OS X, or Windows

### From PyPI (not yet available)

    pip install impy

The package will be available as a pre-compiled binary wheels in the future, but for now you have to compile it from source, see next subsection.

### From source

Installation from source requires a Python installation setup for development, as well as C and Fortran compilers.

To build from source (the **recursive** flag is important to check out submodules):

    git clone --recursive https://github.com/impy-project/impy
    cd impy
    pip install --no-build-isolation --prefer-binary -v -e .

This takes a while. The command `pip install -v -e .` installs the package in editable mode (for developing the Python layer) and with verbose output, so that you can watch the compilation happening. Warnings can be ignored, but watch out for errors.

To run the tests or try the examples, it is convenient use this modified `pip install` instead:

    pip install -v -e .'[test,examples]'

This installs impy and additional optional Python packages to run tests and examples.

#### Known issues

- On OSX
    - You need to install gcc and gfortran with homebrew, and you need to set those compilers explicitly. If you have gcc-12, do `CC=gcc-12 CXX=g++-12 pip install ...`
    - Apple introduced a bug in the Xcode Command Line Tools Version 14 which produces a linker error when compiling C++ code with gcc. Until this is fixed, the workaround is to downgrade to 13.4, use this link https://download.developer.apple.com/Developer_Tools/Command_Line_Tools_for_Xcode_13.4/Command_Line_Tools_for_Xcode_13.4.dmg and turn off automatic updates in the System Settings, because otherwise your Mac will upgrade to 14 again.
- setuptools > 60 does not seem to work. Downgrade with `pip install setuptools<60` if you experience problems.

If you cannot fix the installation with these hints, please look into the subsection below which explains how to install in impy in a verified docker environment. The docker environment has a properly set up environment verified by us, so that the installation is guaranteed to succeed.

### From source in Docker

This guide works on Linux and OSX. You need a running Docker server. Please google how to set up Docker on your machine.

    # download impy
    git clone --recursive https://github.com/impy-project/impy
    cd impy

    # download linux image for x86_64 or see below
    docker pull quay.io/pypa/manylinux2014_x86_64
 
    # For aarch64 or VM on Apple Silicon use the following image and
    # replace the end of the next command accordingly.
    # docker pull quay.io/pypa/manylinux2014_aarch64
    
    # create docker instance and bind impy directory
    docker run --rm -d -it --name impy -v "$(pwd)":/app quay.io/pypa/manylinux2014_x86_64

    # enter your docker instance
    docker exec -it impy /bin/bash

    cd /app

    # select python version, e.g. 3.9, and enter virtual environment
    python3.9 -m venv venv
    source venv/bin/activate

    # install impy and dependencies (prefer binary wheels for deps)
    pip install --prefer-binary -v -e .

You can now use impy inside the docker instance. If you run Linux, you can also make a wheel inside
docker and install it in your host.

    # inside docker
    pip install wheel
    python setup.py bdist_wheel

    # exit docker with ctrl+D
    pip install dist/*.whl

This should allow you to use impy also outside docker. This works only if you use the same Python version inside and outside of docker.

## User interface

There are two ways to interact with the code.

1. As in the example above, via plain python in scripts or Jupyter notebooks. Look at this [example](examples/compare_models.ipynb).

2. Via a HEPMC output that can be piped in Rivet or other tools supporting the format.

## Supported models

- DPMJET-III 3.0.6
- DPMJET-III 19.1
- EPOS-LHC
- PHOJET 1.12-35
- PHOJET 19.1
- PYTHIA 6.4
- PYTHIA 8.3
- QGSJet-01
- QGSJet-II-03
- QGSJet-II-04
- SIBYLL-2.1
- SIBYLL-2.3
- SIBYLL-2.3c
- SIBYLL-2.3d
- SOPHIA (needs update)
- UrQMD 3.4

## Authors:

- Anatoli Fedynitch
- Hans Dembinski
- Anton Prosekin
- Sonia El Hadri
- Keito Watanabe

## LICENSE

The source code of impy is licensed under the [BSD 3-clause license (see LICENSE for detail)](LICENSE). The source codes of the event generators are individually licensed under different conditions (see the COPYING files located in the subdirectories). 
