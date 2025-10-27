# PythiaCascadePython: Pure Python Implementation

This directory contains a pure Python implementation of the PythiaCascade functionality for hadron-nucleus interactions, equivalent to the C++ PythiaCascade class but implemented entirely in Python using chromo's existing Pythia8 integration.

## Files

### Core Implementation
- `pythia_cascade_python.py` - Main implementation of the PythiaCascadePython class

### Examples and Tests
- `test_pythia_cascade_python.py` - Basic functionality tests
- `atmospheric_cascade_example.py` - Atmospheric shower simulation example

## Features

The `PythiaCascadePython` class provides:

1. **Cross Section Calculations**
   - `sigma_setup_hN()` - Calculate hadron-nucleon cross sections
   - `sigma_hA()` - Calculate hadron-nucleus cross sections
   - Nucleus collision multiplicity interpolation

2. **Collision Generation**
   - `next_coll()` - Generate hadron-nucleus collision events
   - Multiple subcollision handling
   - Proper event record management

3. **Physics Implementation**
   - Geometric series probability for multiple collisions
   - Leading particle selection logic
   - Target nucleon (proton/neutron) sampling
   - Energy-momentum conservation

## Key Advantages

- **Pure Python**: No C++ compilation required
- **Transparent**: All logic visible and modifiable
- **Integrated**: Uses chromo's existing event handling
- **Flexible**: Easy to customize physics models

## Usage Example

```python
from chromo.models.pythia_cascade_python import PythiaCascadePython
from chromo.kinematics import FixedTarget
from chromo.constants import GeV

# Set up 1 TeV proton on nitrogen-14
evt_kin = FixedTarget(1000 * GeV, "p", (7, 14))

# Create cascade generator
cascade = PythiaCascadePython(
    evt_kin,
    seed=42,
    max_energy=1e6,  # 1 PeV maximum
    rapid_decays=True,
)

# Generate events
for event in cascade(10):
    print(f"Generated {len(event)} particles")
    # Analyze event...
```

## Physics Implementation Details

### Cross Section Calculation
The implementation uses the same nucleus data tables and interpolation algorithms as the C++ version:
- Tabulated values for A = 1, 2, 4, 9, 12, 14, 16, 27, 40, 56, 63, 84, 107, 129, 197, 208
- Power-law interpolation for intermediate mass numbers
- Two-regime cross section dependence (below/above 20 mb)

### Collision Generation
The `next_coll()` method implements the full collision logic:
1. Initialize event record with incoming hadron
2. Loop over potential nucleon targets with geometric series probability
3. Select leading particle from previous subcollision
4. Choose target nucleon (proton or neutron) 
5. Perform individual hadron-nucleon collision using Pythia8
6. Insert target nucleon and secondary particles
7. Update particle history and relationships

### Event Record Management
- Proper status codes (-181/-182 for target nucleons)
- Mother-daughter relationship tracking
- Energy-momentum bookkeeping
- Optional final state compression

## Limitations

1. **Performance**: Python implementation is slower than C++
2. **Cross Section Accuracy**: Depends on Pythia8's cross section calculations
3. **Simplified Physics**: Some advanced features may need additional implementation

## Testing

Run the test suite:
```bash
python test_pythia_cascade_python.py
```

Run the atmospheric cascade example:
```bash
python atmospheric_cascade_example.py
```

## Technical Notes

### Nucleus Data Tables
The implementation includes the complete nucleus data tables from PythiaCascade.h:
- `_tab_A`: Mass numbers [1, 2, 4, ..., 208]
- `_tab_offset`: Offset parameters for high cross sections
- `_tab_slope`: Slope parameters for high cross sections  
- `_tab_slope_lo`: Slope parameters for low cross sections
- `_tab_border`: Transition point at 20 mb

### Constants
- `_prob_sd = 0.3`: Fraction of single diffractive events
- `_e_kin_min = 0.2`: Minimum kinetic energy threshold (GeV)
- `_mp = 0.938272`: Proton mass (GeV)

### Event Data Structure
The event data is stored as dictionaries with numpy arrays:
```python
{
    "pid": np.array([...]),      # PDG particle IDs
    "status": np.array([...]),   # Status codes
    "charge": np.array([...]),   # Electric charges
    "px", "py", "pz": ...,       # 3-momentum components
    "en": np.array([...]),       # Energies
    "m": np.array([...]),        # Masses
    "vx", "vy", "vz", "vt": ..., # 4-vertex components
    "mothers": [...],            # Mother particle indices
    "daughters": [...],          # Daughter particle indices
}
```

## Future Improvements

1. **Performance Optimization**: Vectorize operations where possible
2. **Enhanced Physics**: Add more sophisticated nuclear models
3. **Cross Section Tables**: Pre-compute and cache cross sections
4. **Decay Handling**: More complete implementation of rapid decays
5. **Validation**: Extensive comparison with C++ PythiaCascade

## References

- PythiaCascade.h - Original C++ implementation
- main483.cc - Atmospheric shower example
- Eur. Phys. J. C82 (2022) 21 (arXiv:2108.03481) - Physics model
