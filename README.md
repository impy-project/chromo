# impy

This package implements are generic user interface to hadronic Interaction Models in PYthon, hence the name *impy*.
The purpose of the package is to liberate particle physicists from Fortran style interfaces to event generators via 'input cards' and
ASCII files. As an open-source project, it is independent of any experimental frameworks or closed-source projects that might restrict your few-author ideas.

To generate some interactions with any of the supported models, it is enough to

```python
from impy.definitions import *
from impy.constants import *
from impy.kinematics import EventKinematics
from impy.common import impy_config, pdata

# Define the parameters of the collisions
event_kinematics = EventKinematics(
    ecm=13 * TeV, p1pdg=2212, p2pdg=2212)

# Create an instance of an event generator by passing
# the model name as a string
generator = make_generator_instance(
    interaction_model_by_tag['SIBYLL2.3C'])

# Initialize it
generator.init_generator(event_kinematics)

# Number of events to generate
nevents = 100

for event in generator.event_generator(event_kinematics, nevents):
    event.filter_final_state_charged()
    # do something with event.p_ids, event.eta, event.en, event.pt, etc.
    # these variables are numpy arrays, that can be histogrammed or counted like
    average_pt = np.mean(event.pt[np.abs(event.p_ids) == 211])
    # This will calculate the average transverse momentum for charged pions

```

## User interface

There are two ways to interact with the code.

1. As in the example above, via plain python in scripts or jupyter notebooks. Look at this [insert link] example.

2. Via a HEPMC output that can be piped in Rivet or other tools supporting the format.



## Supported models

- SIBYLL-2.1
- SIBYLL-2.3
- SIBYLL-2.3c
- EPOS-LHC
- QGSJet-II-03
- QGSJet-II-04
- QGSJet-01
- DPMJET-III 3.0.6
- DPMJET-III 19.1
- PHOJET 1.12-35
- PHOJET 19.1
- UrQMD 3.4

and of course

- PYTHIA 8

## Requirements

See requirements.txt.

## Contributers:

Hans Dembinski
Sonia El Hadri
Keito Watanabe

## LICENSE

[BSD 3-clause license](LICENSE)
