# impy - (hadronic) interaction models in python

This package implements are generic user interface to popular event generators used in cosmic ray and high-energy particle physics. The purpose of the package is to simplify working with simulations of particle interactions without the need to use Fortran style interfaces to event generators, 'ASCII input cards' and files or C++ dependencies.  

Simulate interactions with one of the supported event generators 

```python
import numpy as np
import impy
from impy.constants import TeV

# Define the parameters of the collisions
event_kinematics = impy.kinematics.EventKinematics(ecm=13 * TeV, p1pdg=2212, p2pdg=2212)
# Create an instance of an event generator
generator = impy.models.Sibyll23d(event_kinematics)

nevents = 0
average_pt = 0

# Generate 10000 events
for event in generator(10000):
    # Filter events
    event.filter_final().filter_charged()
    # do something with event.p_ids, event.eta, event.en, event.pt, etc.
    # these variables are numpy arrays, that can be histogrammed or counted like
    pt = event.pt[np.abs(event.p_ids) == 211]
    # The list could be empty
    if len(pt) > 0:
        nevents += 1
        average_pt += np.mean(pt)

average_pt = 1 / float(nevents) * average_pt
print("Average pT for charged pions {0:4.3f}".format(average_pt))
```

## Installation

The package is (will be) available including pre-compiled binaries. The installation in that case simplifies to (*this does not work yet use installation from source*):

    pip install impy

To build from source (the **recursive** flag is important to checkout the sub-modules):

    git clone --recursive https://github.com/afedynitch/impy
    cd impy
    pip install -e .
    make -j<insert number of CPU cores>

For now, you need to call make by hand, but this will be automated. The command `pip install -e .` installs the package in editable mode (for developers).

Because of the architectural transition and there are many issues on mac, building from source may be a bit complicated. Using brew gcc and python it is possible to build the code by:

    CC=gcc-10 CXX=gcc-10 FC=gfortran-10 PYTHON_EXE=/usr/local/opt/python@3.8/bin/python3 make -jXXX

Replace `gcc-10` by your version in brew. The official Mac Python is currently broken due to th transition to Apple Silicon, but it is possible to build with a bit of hacking. But currently
I don't use a Mac and cannot debug it. 
 
## Requirements

- Python 2.7 - 3.9
- Linux, Mac OS X, or Windows
- pip
- particletools
- numpy
- scipy
- pyyaml
- pyhepmc

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
- Sonia El Hadri
- Keito Watanabe

## LICENSE

The source code of impy is licensed under the [BSD 3-clause license (see LICENSE for detail)](LICENSE). The source codes of the event generators are individually licensed under different conditions (see the COPYING files located in the subdirectories). 
