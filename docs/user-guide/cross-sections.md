# Cross Sections

Every generator can compute cross sections for the current kinematics.

## Basic Usage

```python
xs = generator.cross_section()
print(f"Total:     {xs.total:.2f} mb")
print(f"Inelastic: {xs.inelastic:.2f} mb")
print(f"Elastic:   {xs.elastic:.2f} mb")
```

## CrossSectionData Fields

All values are in **millibarn** (mb). Fields not provided by a specific model are set to `NaN`.

| Field | Description |
|-------|-------------|
| `total` | Total cross section (elastic + inelastic) |
| `inelastic` | Inelastic cross section (includes diffractive) |
| `elastic` | Elastic scattering cross section |
| `prod` | Production cross section (total - elastic) |
| `quasielastic` | Quasielastic cross section (includes elastic) |
| `coherent` | Coherent (elastic w.r.t. projectile) cross section |
| `diffractive_xb` | Single diffractive: target intact |
| `diffractive_ax` | Single diffractive: projectile intact |
| `diffractive_xx` | Double diffractive |
| `diffractive_axb` | Central diffractive |
| `diffractive_sum` | Sum of diffractive components |
| `b_elastic` | Elastic slope in mb/GeV^2 |
| `emd` | Electromagnetic dissociation cross section |

!!! note
    The level of detail varies by model. Pythia 8 fills most fields, while simpler models may only provide `total` and `inelastic`. Check for `NaN` before using a field.

## Cross Sections at Different Kinematics

You can query cross sections without changing the generator's current state:

```python
from chromo.kinematics import CenterOfMass
from chromo.constants import GeV

# Query at a different energy (does not change generator state)
xs_low = generator.cross_section(CenterOfMass(50 * GeV, "proton", "proton"))
xs_high = generator.cross_section(CenterOfMass(1000 * GeV, "proton", "proton"))
```

## Composite Target Cross Sections

For `CompositeTarget` kinematics, the cross section is the weighted average over the components:

```python
from chromo.kinematics import CompositeTarget

air = CompositeTarget([("N14", 0.7843), ("O16", 0.2105), ("Ar40", 0.0052)])
kin = CenterOfMass(100 * GeV, "proton", air)
xs = generator.cross_section(kin)
```
