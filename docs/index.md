# chromo

**Cosmic Ray and Hadronic Interaction Monte Carlo Frontend**

chromo provides a simple, unified Python interface to popular hadronic event generators used in cosmic-ray and high-energy particle physics. It removes the need for Fortran-style interfaces, ASCII input cards, and complex C++ dependencies.

## Quick Example

```python
import chromo

kinematics = chromo.kinematics.CenterOfMass(
    13 * chromo.constants.TeV, "proton", "proton"
)
generator = chromo.models.Sibyll23d(kinematics)

for event in generator(100):
    event = event.final_state_charged()
    print(f"Event {event.nevent}: {len(event)} charged particles")
```

## Supported Generators

chromo wraps **9 model families** with **26 model classes**, covering hadron-hadron, hadron-nucleus, nucleus-nucleus, photon-hadron, and electron-positron collisions:

DPMJET-III, EPOS-LHC, FLUKA, PHOJET, PYTHIA 6, PYTHIA 8 (including Cascade and Angantyr), QGSJet, SIBYLL, SOPHIA, and UrQMD.

See the [Model Overview](models/overview.md) for a full capability table.

## Getting Started

- **[Installation](getting-started/installation.md)** -- install from PyPI or build from source
- **[Your First Simulation](getting-started/first-simulation.md)** -- step-by-step walkthrough
- **[User Guide](user-guide/kinematics.md)** -- kinematics, events, cross sections, and more

## Citation

If you use chromo in your research, please cite:

> A. Fedynitch, H. Dembinski and A. Prosekin, *Chromo: A high-performance python interface to hadronic event generators for collider and cosmic-ray simulations*, [Comput.Phys.Commun. 321 (2026) 110031](https://doi.org/10.1016/j.cpc.2026.110031), [arXiv:2507.21856](https://arxiv.org/abs/2507.21856)
