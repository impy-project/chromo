# Pythia8Cascade Implementation

This document describes the new `Pythia8Cascade` class that extends the chromo library to support hadron-nucleus interactions using the PythiaCascade functionality from Pythia8.

## Overview

The `Pythia8Cascade` class is a derived class from the base `Pythia8` class that implements the PythiaCascade wrapper for simulating hadron-nucleus collisions. This is particularly useful for:

- Cosmic ray simulations
- High-energy hadron-nucleus interactions
- Fixed-target experiments with nuclear targets
- Cascade simulations in dense media

## Key Features

1. **Multiple Energy Support**: Can handle interactions at various energies up to a specified maximum
2. **Nuclear Target Support**: Supports a wide range of nuclear targets from deuterium to lead
3. **MPI Table Acceleration**: Can use pre-computed MPI tables (like main424.cmnd) for faster initialization
4. **Fixed Target Framework**: Operates in fixed-target mode for realistic hadron-nucleus collisions

## Usage

### Basic Initialization

```python
from chromo.models.pythia8_cascade import Pythia8Cascade
from chromo.kinematics import FixedTarget

# Create kinematics for proton-carbon collision at 1 TeV
kin = FixedTarget(1000, "p", "C12")  # 1 TeV proton on carbon

# Initialize cascade with default settings
cascade = Pythia8Cascade(kin, seed=42)
```

### Advanced Configuration

```python
# Initialize with custom settings
cascade = Pythia8Cascade(
    kin, 
    seed=42,
    max_energy=1e4,  # Maximum energy in GeV for cascade initialization
    use_mpi_tables=True,  # Use pre-computed MPI tables for acceleration
    mpi_file="pythiaCascade.mpi"  # MPI table file name
)
```

### Supported Target Nuclei

The class supports the following nuclear targets (among others):
- Light nuclei: He4, Li6, Be9, C12, N14, O16
- Medium nuclei: Al27, Ar40, Fe56, Cu63, Kr84
- Heavy nuclei: Ag107, Xe129, Au197, Pb208

### Event Generation

```python
# Generate cross section information
xs = cascade.cross_section()
print(f"Total cross section: {xs.total:.2f} mb")
print(f"Inelastic cross section: {xs.inelastic:.2f} mb")

# Generate events
for event in cascade(100):  # Generate 100 events
    print(f"Event with {len(event.pid)} particles")
    
    # Access final state particles
    final_particles = event.pid[event.status == 1]
    print(f"Final state particles: {len(final_particles)}")
```

### Multiple Target Types

```python
# Different nuclear targets
targets = ["C12", "Fe56", "Pb208"]
energies = [100, 1000, 10000]  # GeV

for target in targets:
    for energy in energies:
        kin = FixedTarget(energy, "p", target)
        cascade = Pythia8Cascade(kin, max_energy=2*energy)
        
        # Calculate cross section
        xs = cascade.cross_section()
        print(f"{energy} GeV p-{target}: σ_total = {xs.total:.2f} mb")
```

## Implementation Details

### PythiaCascade Integration

The class wraps the C++ `PythiaCascade` functionality through pybind11 bindings:

1. **Initialization**: Uses `PythiaCascade::init()` with appropriate energy and MPI settings
2. **Cross Section Calculation**: Implements `sigmaSetuphN()` and `sigmahA()` for hadron-nucleus cross sections
3. **Event Generation**: Uses `nextColl()` for hadron-nucleus collisions

### C++ Binding Extensions

The pybind11 interface (`_pythia8.cpp`) has been extended to include:
- `PythiaCascade` class binding
- Support for 4-momentum and vertex arrays
- Event generation methods (`nextColl`, `nextDecay`)
- Cross section calculation methods

### Acceleration Features

1. **MPI Tables**: Can use pre-computed multiparton interaction tables from `main424.cmnd`
2. **Energy Optimization**: Initializes with maximum expected energy to avoid re-initialization
3. **Efficient Cross Sections**: Caches hadron-nucleon cross sections for multiple A values

## Differences from Base Pythia8 Class

| Feature | Pythia8 | Pythia8Cascade |
|---------|---------|----------------|
| Frame | Center of Mass | Fixed Target |
| Target Support | Limited nuclei | Extended nuclear targets |
| Physics Focus | Elementary collisions | Hadron-nucleus cascades |
| Energy Range | Broad | Optimized for high energy |
| Initialization | Per-event | Energy-range optimized |

## Technical Notes

### Memory Management
- Events are stored temporarily in the generator for access by the event wrapper
- C++ objects are managed through pybind11 with appropriate reference policies

### Performance Considerations
- Initialize with the highest energy expected in your simulation
- Use MPI tables when available for faster startup
- The cascade approach is most efficient for high-energy interactions

### Limitations
- Currently focuses on inelastic processes (elastic may be added later)
- Charge calculation is simplified (full implementation pending)
- Some advanced PythiaCascade features not yet exposed

## Example Applications

### Cosmic Ray Simulation
```python
# High-energy cosmic ray proton hitting atmosphere
kin = FixedTarget(1e6, "p", "N14")  # 1 PeV proton on nitrogen
cascade = Pythia8Cascade(kin, max_energy=1e7)

for event in cascade(1000):
    # Analyze atmospheric shower development
    analyze_shower_particles(event)
```

### Fixed Target Experiment
```python
# Beam-target experiment
beam_energy = 158  # GeV (SPS energy)
kin = FixedTarget(beam_energy, "p", "Pb208")
cascade = Pythia8Cascade(kin, max_energy=200)

# Study particle production
for event in cascade(10000):
    analyze_particle_production(event)
```

This implementation provides a powerful tool for hadron-nucleus cascade simulations while maintaining compatibility with the existing chromo framework.
