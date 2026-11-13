# DPMJET & PHOJET

DPMJET (Dual Parton Model with JETs) and PHOJET (PHOton JETs) are closely related model families.
PHOJET handles photon-induced interactions and is embedded as a component of the DPMJET-III
codebase: in versions 19.x, PHOJET is the photon-interaction engine inside DPMJET-III. In the
3.0.7 generation, they share the same underlying two-component Dual Parton Model framework but
are packaged as separate model classes in chromo.

---

## Version matrix

| DPMJET class | PHOJET class | Generation | Nuclear A | Extended projectiles | PHOJET targets |
|---|---|---|---|---|---|
| `DpmjetIII307` | `Phojet112` | Legacy | up to 280 | No | p, n |
| `DpmjetIII191` | `Phojet191` | Modern | up to 280 | Yes | p, n, γ (γγ, γN) |
| `DpmjetIII193` | `Phojet193` | Latest | up to 280 | Yes | p, n, γ (γγ, γN) |

**Recommendation:** use 19.3 for new projects. Use 3.0.7 only for reproducing older results
or comparisons to published calculations based on that version.

---

## DPMJET — nuclear collisions

### Supported systems

All DPMJET versions support nucleus–nucleus collisions up to **A = 280** on both beam and target,
covering the full range from proton–proton to Pb+Pb and Au+Au. The minimum center-of-mass energy
is `_ecm_min = 1 GeV`.

Targets are specified as nuclei using the `(Z, A)` tuple or standard name (`"Pb208"`, `"Au197"`,
etc.); the `Nuclei()` helper class provides convenience constructors.

### Extended projectile set (19.x only)

DPMJET 19.1 and 19.3 add hadronic projectiles beyond the standard p, n, π±, K± set:

- K-mesons: K⁰_S, K⁰_L
- D-mesons: D⁰, D±, D_s
- Lambda, Sigma, Xi, Omega baryons and their charge conjugates

These are available in both `DpmjetIII191` and `DpmjetIII193`.

### Heavy-ion example (Pb+Pb at LHC energy)

```python
from chromo.models import DpmjetIII193
from chromo.kinematics import CenterOfMass

gen = DpmjetIII193(CenterOfMass(5020, "Pb208", "Pb208"), seed=1)

for event in gen(10):
    fs = event.final_state_charged()
    print(f"Nch = {len(fs)}")
```

### Proton–nucleus example

```python
from chromo.models import DpmjetIII193
from chromo.kinematics import FixedTarget

# 100 TeV proton on iron-56, fixed target
gen = DpmjetIII193(FixedTarget(1e5, "p", "Fe56"), seed=42)

for event in gen(50):
    fs = event.final_state()
    print(f"N = {len(fs)}, E_sum = {fs.en.sum():.1f} GeV")
```

### Cross sections

```python
xs = gen.cross_section()
print(f"σ_total   = {xs.total:.1f} mb")
print(f"σ_inel    = {xs.inelastic:.1f} mb")
print(f"σ_elastic = {xs.elastic:.1f} mb")
```

---

## PHOJET — photon-induced collisions

PHOJET implements the two-component Dual Parton Model for photon interactions. The physics scope
and supported targets differ by version:

| Version | Supported systems |
|---------|------------------|
| `Phojet112` | γ+p, γ+n |
| `Phojet191` | γ+p, γ+n, γ+γ |
| `Phojet193` | γ+p, γ+n, γ+γ |

### Photon–proton example

```python
from chromo.models import Phojet193
from chromo.kinematics import CenterOfMass

gen = Phojet193(CenterOfMass(200, "gamma", "p"), seed=7)

for event in gen(100):
    fs = event.final_state()
    print(f"N = {len(fs)}")
```

### Gamma–gamma example (19.x only)

```python
from chromo.models import Phojet193
from chromo.kinematics import CenterOfMass

gen = Phojet193(CenterOfMass(500, "gamma", "gamma"), seed=3)

for event in gen(50):
    fs = event.final_state_charged()
    print(f"Nch = {len(fs)}")
```

---

## Gotchas and limitations

!!! danger "Single instantiation per process"
    Each DPMJET/PHOJET version uses Fortran COMMON blocks and can only be initialized once per
    Python process. Different versions (e.g., `DpmjetIII191` and `DpmjetIII193`) cannot coexist
    in the same process.

!!! note "DPMJET energy floor"
    `_ecm_min = 1 GeV` applies to all DPMJET variants. Below this threshold, kinematics
    validation will reject the configuration.

!!! note "PHOJET 1.12 target restriction"
    `Phojet112` does not support photon beams on photon targets. Use `Phojet191` or `Phojet193`
    for γγ collisions.

!!! note "Nuclear targets in PHOJET"
    PHOJET treats only free nucleon (proton, neutron) and photon targets. It does not support
    nuclear targets (A > 1). For photon–nucleus interactions, use a model with nuclear support.

!!! info "DPMJET initialization time"
    DPMJET 19.x initializes over 700 Fortran source files and may take several seconds to start.
    Plan for a noticeable cold-start delay when first creating a generator instance.
