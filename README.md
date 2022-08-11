# impy - (hadronic) interaction models in python

This package implements are generic user interface to popular event generators used in cosmic ray and high-energy particle physics. The purpose of the package is to simplify working with simulations of particle interactions without the need to use Fortran style interfaces to event generators, 'ASCII input cards' and files or C++ dependencies.  

Simulate interactions with one of the supported event generators 

```python
from impy.definitions import *
from impy.constants import *
from impy.kinematics import EventKinematics
from impy import impy_config

# Define the parameters of the collisions
event_kinematics = EventKinematics(
    ecm=13 * TeV, p1pdg=2212, p2pdg=2212)

# Create an instance of an event generator by passing
# the model name as a string
generator = make_generator_instance(
    interaction_model_by_tag['SIBYLL2.3D'])

# Initialize it
generator.init_generator(event_kinematics)

# Number of events to generate
nevents = 100

for event in generator.event_generator(event_kinematics, nevents):
    event.filter_final_state_charged()
    # do something with event.p_ids, event.eta, event.en, event.pt, etc.
    # these variables are numpy arrays, that can be histogrammed or counted like
    average_pt += 1/float(nevents)*np.mean(event.pt[np.abs(event.p_ids) == 211])

print('Average pT for charged pions {0:4.3f}'.format(average_pt))
```

## Installation

## Supported platforms

- Python 3.6+
- Linux, Mac OS X, or Windows

### Without docker

If you have trouble with this installation guide, look into the subsection which explains how to install in impy in a fixed docker environment.

The package is (will be) available including pre-compiled binaries. The installation in that case simplifies to (*this does not work yet use installation from source*):

    pip install impy

To build from source (the **recursive** flag is important to checkout the sub-modules):

    git clone --recursive https://github.com/impy-project/impy
    cd impy
    pip install -e .
    make -j<insert number of CPU cores>

For now, you need to call make by hand, but this will be automated. The command `pip install -e .` installs the package in editable mode (for developers).

Because of the architectural transition and there are many issues on mac, building from source may be a bit complicated. Using brew gcc and python it is possible to build the code by:

    CC=gcc-10 CXX=gcc-10 FC=gfortran-10 PYTHON_EXE=/usr/local/opt/python@3.8/bin/python3 make -jXXX

Replace `gcc-10` by your version in brew. The official Mac Python is currently broken due to the transition to Apple Silicon, but it is possible to build with a bit of hacking. But currently
I don't use a Mac and cannot debug it. 
 
### With docker

This guide works on Linux and OSX. You need a running docker server. Please google how to set up docker on your machine.

    # download impy
    git clone --recursive https://github.com/impy-project/impy
    cd impy

    # download linux image for x86_64 or see below
    docker pull quay.io/pypa/manylinux2014_x86_64
 
    # For aarch64 or VM on Apple Silicon use the following image and
    # replace the end of the next command accordingly.
    # docker pull quay.io/pypa/manylinux2014_aarch64
    
    # create docker instance and bind impy directory
    docker run -d -it --name impy -v "$(pwd)":/app quay.io/pypa/manylinux2014_x86_64

    # enter your docker instance
    docker exec -it impy /bin/bash

    cd /app

    # select python version, e.g. 3.8, and enter virtual environment
    python3.8 -m venv venv
    source venv/bin/activate

    # install impy and dependencies (prefer binary wheels for deps)
    pip install --prefer-binary -e .

    # compile the FORTRAN interface (this will be automated in the future)
    make -j<insert number of CPU cores>

You can now use impy inside the docker instance.

## User interface

There are two ways to interact with the code.

1. As in the example above, via plain python in scripts or jupyter notebooks. Look at this [example](examples/compare_two_models.ipynb).

2. Via a HEPMC output that can be piped in Rivet or other tools supporting the format.

## Supported models

- DPMJET-III 3.0.6
- DPMJET-III 19.1
- EPOS-LHC
- PHOJET 1.12-35
- PHOJET 19.1
- PYTHIA 6
- PYTHIA 8 (not yet bundled)
- QGSJet-01
- QGSJet-II-03
- QGSJet-II-04
- SIBYLL-2.1
- SIBYLL-2.3
- SIBYLL-2.3c
- SIBYLL-2.3d
- SOPHIA (needs update)
- DPMJET-II (also needs update but model deprecated)
- UrQMD 3.4


## Authors:

- Anatoli Fedynitch
- Hans Dembinski
- Anton Prosekin
- Sonia El Hadri
- Keito Watanabe

## LICENSE

The source code of impy is licensed under the [BSD 3-clause license (see LICENSE for detail)](LICENSE). The source codes of the event generators are individually licensed under different conditions (see the COPYING files located in the subdirectories). 
