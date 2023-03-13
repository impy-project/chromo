# ![](doc/chromo.svg)<br> Cosmic ray and HadROnic interactiOn MOnte-carlo frontend

This package provides a simple and generic user interface to popular event generators used in cosmic ray and high-energy particle physics. By removing the need for complicated Fortran-style interfaces, ASCII input cards, and C++ dependencies, the package simplifies the simulation of particle interactions, making it easier and faster for a wider audience to access.

## Usage

### Python user interface

To simulate interactions with one of the supported event generators, import the package and define the parameters of the collision. Then, create an instance of an event generator, and generate events.

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

Further examples, such as [this](examples/compare_models.ipynb) can be found in the examples folder.

### Command line user interface (CLI) 

CLI via a HEPMC output that can be piped in Rivet or other tools supporting the format.

## Output formats

- plain Python `event` objects
- HEPMC (optionally gzip compressed)
- ROOT (via uproot)
- event views via SVG  

## Supported models

Please note that `chromo` only provides a user interface for the following models, and does not contain any particle physics models itself. When using any of these models in public-facing work, it is important to properly cite the original model reference by following the links below. Additionally, if you find `chromo` useful in your work, we would appreciate an acknowledgement, footnote, or link to `chromo`.

| Interaction model                                         | Supported proj/targ       | Comment                         | 
|------------------------------------------------------------|---------------------------|--------------------------------|
| [DPMJET-III 3.0.6](https://inspirehep.net/literature/538940) & [PHOJET 1.12-35](https://inspirehep.net/literature/373339)      | *hN, gg, gN, hA, gA, AA*  | |
| [DPMJET-III & PHOJET 19.1 and 19.3](https://inspirehep.net/literature/1503512) [(repo on GitHub)](https://github.com/DPMJET/DPMJET) |  *hN, gg, gN, hA, gA, AA* | |
| [EPOS-LHC](https://inspirehep.net/literature/1236629)     | *hN, hA, AA*              | |
| [PYTHIA 6.4](https://inspirehep.net/literature/712925)    | *hN, ee, gg, gN*          | |
| [PYTHIA 8.3](https://inspirehep.net/literature/2056998) (https://pythia.org/) | *hN, ee, gg, gN* & *hA, AA* (Argantyr) | unavailable on Windows |
| [QGSJet-01](https://inspirehep.net/literature/460408)     | *hN, hA, AA*              | |
| [QGSJet-II-03](https://inspirehep.net/literature/667881)  | *hN, hA, AA*              | |
| [QGSJet-II-04](https://inspirehep.net/literature/872658)  | *hN, hA, AA*              | |
| [SIBYLL-2.1](https://inspirehep.net/literature/823839)    | *hN, hA (A<=20)*          | |
| [SIBYLL-2.3d](https://inspirehep.net/literature/1768983)  | *hN, hA (A<=20)*          | incl. legacy versions -2.3/-2.3c |
| [SOPHIA 2.0](https://inspirehep.net/literature/497602)    | *gN*                      | |
| [UrQMD 3.4](https://inspirehep.net/literature/468266) [+ second citation](https://inspirehep.net/literature/507334)    |  hN, hA, AA* | unavailable on Windows |


*h* = hadron, *A* = nucleus, *g* = gamma, *e* = electron/positron

## Installation via PyPI

### Supported platforms

- Python 3.8+
- Linux, Mac OS X (x86 and M1/M2), Windows

The recommended way to install `chromo` is by using the pre-compiled binary wheel, which is available for most common architectures and Python versions

    pip install chromo

Advanced and developer installation instructions can be found [here](doc/dev_docs.md).


## Authors

- Anatoli Fedynitch
- Hans Dembinski
- Anton Prosekin
- Sonia El Hadri
- Keito Watanabe

## LICENSE

The source code of chromo is licensed under the [BSD 3-clause license (see LICENSE for detail)](LICENSE). The source codes of the event generators are individually licensed under different conditions (see the COPYING files located in the subdirectories).
