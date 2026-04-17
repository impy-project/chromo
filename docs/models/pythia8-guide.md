# Pythia 8 Guide

chromo exposes three Pythia 8 model classes with distinct physics scopes. Choosing the right one
depends on your collision system and energy range.

!!! warning "Windows not supported"
    All three Pythia 8 classes are **unavailable on Windows**. They are conditionally excluded
    from the build when `host_machine.system() == 'windows'`.

---

## Decision tree

```
What is your collision system?
│
├── hN, ee, γγ, γN — and no nuclear targets needed
│   └── Use Pythia8
│       (also supports extended projectiles: K+, D+, B+, Lambda, Xi, ...)
│
├── h+A, single-collision approximation
│   └── Use Pythia8Cascade
│       (nuclear projectiles decomposed to Z p + (A-Z) n; target A > 1)
│
└── hA or AA, full Glauber geometry
    └── Use Pythia8Angantyr
        (precomputed tables 20 GeV–20 PeV CMS; target must be nucleus A > 1)
```

---

## `Pythia8` — standard hadron/lepton collisions

**Scope:** proton–proton, proton–antiproton, e⁺e⁻, γγ, γN. Extended projectile set covers
strange, charm, and bottom hadrons (K±, K⁰, D±, D⁰, B±, B⁰, Λ, Σ, Ξ, Ω and their
antiparticles). **No nuclear targets.**

```python
from chromo.models import Pythia8
from chromo.kinematics import CenterOfMass

gen = Pythia8(CenterOfMass(13000, "p", "p"), seed=42)

for event in gen(100):
    fs = event.final_state_charged()
    print(f"Nch = {len(fs)}")
```

### Cross sections

`cross_section()` returns the full Pythia 8 cross-section breakdown including total, elastic,
inelastic, and diffractive components.

```python
xs = gen.cross_section()
print(f"σ_inel = {xs.inelastic:.2f} mb")
```

---

## `Pythia8Cascade` — single-collision h+A

**Scope:** hadronic projectile on a nuclear target (A > 1). Implements the PythiaCascade plugin,
which decomposes nuclear projectiles into Z protons + (A−Z) neutrons and collides them one at a
time. Supports `CompositeTarget` via automatic event splitting over target components.

**Key physics convention:** `slowDecays=True` — follows the cosmic-ray convention and decays
particles that are normally stable in collider contexts (Σ⁰, Ξ⁰, etc.).

```python
from chromo.models import Pythia8Cascade
from chromo.kinematics import FixedTarget

# 1 PeV proton on oxygen-16 in fixed-target mode
gen = Pythia8Cascade(FixedTarget(1e6, "p", (8, 16)), seed=42)

for event in gen(50):
    fs = event.final_state()
    print(f"N_final = {len(fs)}, pids = {fs.pid[:5]}")
```

### CompositeTarget

```python
from chromo.kinematics import CompositeTarget

air = CompositeTarget([("N14", 0.78), ("O16", 0.21), ("Ar40", 0.01)])
gen = Pythia8Cascade(FixedTarget(1e5, "p", air), seed=1)
```

Events are generated for each component proportionally to its weight; the composite event record
is assembled transparently.

### Cross sections

`cross_section()` returns a fast parametric estimate based on the nCollAvg formula — suitable for
rate calculations but not as precise as a full Glauber MC.

```python
xs = gen.cross_section()
print(f"σ_inel(p+O16) ≈ {xs.inelastic:.1f} mb")
```

---

## `Pythia8Angantyr` — Glauber heavy-ion model

**Scope:** hA and AA collisions with full Glauber geometry. Initialization requires precomputed
tables covering **20 GeV to 20 PeV center-of-mass energy** (`_ecm_min = 20 GeV`). Targets must
be nuclei — free proton or neutron targets are **not supported**.

Live target switching is available via `setBeamIDs`, so a single generator instance can be
reused across multiple target species without full re-initialization. Supports `CompositeTarget`.

```python
from chromo.models import Pythia8Angantyr
from chromo.kinematics import CenterOfMass

gen = Pythia8Angantyr(CenterOfMass(5020, "p", "Pb208"), seed=7)

for event in gen(20):
    fs = event.final_state()
    print(f"N_final = {len(fs)}")
```

### Cross sections — two methods

| Method | Speed | Precision |
|--------|-------|-----------|
| `cross_section()` | Fast | Parametric (nCollAvg formula) |
| `glauber_cross_section(n_trials)` | Slow | Full GlauberOnly MC |

```python
# Fast estimate (same formula as Cascade)
xs_fast = gen.cross_section()

# Precise Angantyr Glauber MC (runs n_trials events)
xs_glauber = gen.glauber_cross_section(n_trials=10000)
print(f"σ_inel(p+Pb) = {xs_glauber.inelastic:.1f} mb")
```

### CompositeTarget with Angantyr

```python
from chromo.kinematics import CompositeTarget

air = CompositeTarget([("N14", 0.78), ("O16", 0.21), ("Ar40", 0.01)])
gen = Pythia8Angantyr(CenterOfMass(5000, "p", air), seed=3)
```

---

## Gotchas and limitations

!!! danger "Single instantiation per process"
    All three Pythia 8 classes share the same compiled Pythia 8 library, which uses global state.
    **Only one instance of any Pythia 8 class can be created per Python process.** Attempting to
    create a second instance (even of a different class) will raise an error. Use separate
    subprocesses or `multiprocessing` if you need multiple configurations.

!!! note "Angantyr energy floor"
    `Pythia8Angantyr` requires `√s ≥ 20 GeV` (the `_ecm_min` limit set by the precomputed
    initialization tables). For lower energies or nuclear projectiles in single-collision
    approximation, use `Pythia8Cascade`.

!!! note "Angantyr target restriction"
    `Pythia8Angantyr` does not support free proton (`p`) or neutron (`n`) as targets. The target
    must be a nucleus with A > 1. For p+p or p+n, use `Pythia8`.

!!! info "Cascade slow-decay convention"
    `Pythia8Cascade` runs with `slowDecays=True`, meaning particles that are treated as stable
    at colliders (Σ⁰, Ξ⁰, etc.) are decayed. This matches the cosmic-ray air-shower convention
    used by CORSIKA and similar frameworks. If you compare multiplicity to collider data,
    account for this difference.

!!! info "Cascade RNG reproducibility"
    `Pythia8Cascade` contains two internal Pythia instances (`pythiaMain` and `pythiaColl`).
    RNG state save/restore (`getRndmState` / `setRndmState`) operates on both, ensuring fully
    reproducible event sequences when the seed is fixed.
