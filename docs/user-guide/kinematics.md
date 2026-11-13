# Kinematics

chromo provides a unified way to specify collision kinematics regardless of which event generator you use. All generators accept the same kinematics objects.

## Energy Units

Import energy units from `chromo.constants`:

```python
from chromo.constants import GeV, TeV, PeV, EeV, MeV
```

These are conversion factors. `13 * TeV` gives 13 TeV in chromo's internal units (MeV).

## Center-of-Mass Frame

Use `CenterOfMass` when you specify the center-of-mass energy and want events output in the CMS frame:

```python
from chromo.kinematics import CenterOfMass
from chromo.constants import TeV

kin = CenterOfMass(13 * TeV, "proton", "proton")

print(f"sqrt(s) = {kin.ecm:.1f} MeV")  # internal units are MeV
print(f"E_lab   = {kin.elab:.1f} MeV")
```

## Fixed Target Frame

Use `FixedTarget` when you specify the projectile energy in the laboratory frame (target at rest):

```python
from chromo.kinematics import FixedTarget
from chromo.constants import PeV

kin = FixedTarget(1 * PeV, "proton", "N14")

print(f"E_lab   = {kin.elab:.1f} MeV")
print(f"sqrt(s) = {kin.ecm:.1f} MeV")
```

For nuclear projectiles or targets, the energy is **per nucleon**.

## Specifying Particles

Particles can be given as:

- **Name strings:** `"proton"`, `"neutron"`, `"pi+"`, `"pi-"`, `"K+"`, `"gamma"`
- **PDG IDs:** `2212` (proton), `211` (pi+), `22` (gamma)
- **Nucleus strings:** `"N14"`, `"O16"`, `"Fe56"`, `"Pb208"`
- **Nucleus (A, Z) tuples:** `(14, 7)` for nitrogen

```python
# These are all equivalent for a proton
CenterOfMass(100 * GeV, "proton", "proton")
CenterOfMass(100 * GeV, 2212, 2212)
CenterOfMass(100 * GeV, "p", "p")
```

## Composite Targets

For mixed materials like air, use `CompositeTarget`:

```python
from chromo.kinematics import CompositeTarget, CenterOfMass
from chromo.constants import TeV

# Air is ~78% nitrogen, 21% oxygen, 1% argon (by volume)
air = CompositeTarget(
    [("N14", 0.7843), ("O16", 0.2105), ("Ar40", 0.0052)]
)

kin = CenterOfMass(100 * TeV, "proton", air)
```

When generating events with a `CompositeTarget`, chromo automatically splits the requested number of events across the components according to their fractions and runs the generator for each component.

## Changing Kinematics

You can change the kinematics of an existing generator by assigning to the `kinematics` property:

```python
generator = chromo.models.Sibyll23d(
    CenterOfMass(100 * GeV, "proton", "proton")
)

# Change to a different energy
generator.kinematics = CenterOfMass(1 * TeV, "proton", "proton")
```

The generator will validate that the new kinematics are within its supported range.
