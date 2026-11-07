# Your First Simulation

This walkthrough shows how to set up a collision, run an event generator, and access particle data.

## 1. Define the Collision

Every simulation starts by specifying the collision kinematics: the energy, projectile, and target.

```python
import chromo

# 13 TeV proton-proton in the center-of-mass frame
kinematics = chromo.kinematics.CenterOfMass(
    13 * chromo.constants.TeV, "proton", "proton"
)
```

`CenterOfMass` specifies that the energy is the center-of-mass energy and that events will be output in the CMS frame. You can also use `FixedTarget` for laboratory-frame kinematics:

```python
# 1 PeV proton on nitrogen in the lab frame
kinematics = chromo.kinematics.FixedTarget(
    1 * chromo.constants.PeV, "proton", "N14"
)
```

Particles can be specified by name (`"proton"`, `"pi+"`, `"N14"`) or by PDG ID.

## 2. Create a Generator

Pass the kinematics to a model class:

```python
generator = chromo.models.Sibyll23d(kinematics)
```

!!! warning "Single instantiation"
    Most generators use Fortran global state and can only be instantiated **once per Python process**. Creating a second instance of the same model will raise an error. See [Running a Generator](../user-guide/running-a-generator.md) for how to work with multiple models.

## 3. Generate Events

The generator is a callable that returns an iterator:

```python
for event in generator(1000):
    # event is an EventData object
    # event.pid, event.pt, event.eta, etc. are numpy arrays
    pass
```

## 4. Filter and Analyze

Use `final_state()` or `final_state_charged()` to select relevant particles:

```python
import numpy as np

for event in generator(1000):
    fs = event.final_state_charged()

    # Transverse momentum of charged pions
    is_pion = np.abs(fs.pid) == 211
    if np.any(is_pion):
        mean_pt = np.mean(fs.pt[is_pion])
        print(f"Event {fs.nevent}: mean pion pT = {mean_pt:.3f} GeV")
```

## 5. Get Cross Sections

```python
xs = generator.cross_section()
print(f"Total:     {xs.total:.2f} mb")
print(f"Inelastic: {xs.inelastic:.2f} mb")
print(f"Elastic:   {xs.elastic:.2f} mb")
```

## What's Next?

- [Kinematics](../user-guide/kinematics.md) -- energy specifications, frame conversions, composite targets
- [Working with Events](../user-guide/working-with-events.md) -- all available particle properties
- [Model Overview](../models/overview.md) -- which generator to use for your physics case
