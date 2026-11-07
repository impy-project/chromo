# FLUKA

FLUKA 2025.1 is an optional, license-restricted backend integrated into chromo via
f2py-exported Fortran helpers linked against the official FLUKA archives
(`libflukahp.a`, `libdpmmvax.a`, `librqmdmvax.a`, `librqmd.a`, `libdpmjet*.a`).
It supports **hN, hA, AA** (cross sections only for nuclear projectiles),
**photohadronic**, **photonuclear**, and **EMD** interactions, with nuclear remnants
included in HEPEVT output.

!!! warning "License-restricted — not in public wheels"
    FLUKA requires a license from [fluka.cern](https://fluka.cern). It is not included
    in PyPI distribution or public CI wheels. You must build chromo from source with
    `$FLUPRO` set. See [Building with FLUKA](../getting-started/building-with-fluka.md).

---

## Overview

FLUKA dispatches to different hadronic models depending on the lab-frame kinetic energy
per nucleon:

| Regime | System | Low-energy model | High-energy model | Default transition |
|--------|--------|--------------------|-------------------|-------------------|
| **hA** | hadron + nucleus | PEANUT | DPMJET-3 | 20 TeV/nucleon (lab kinetic) |
| **AA** | nucleus + nucleus | rQMD (above ~125 MeV/n) | DPMJET-3 | 12.5 GeV/n ± 2.0 GeV/n |
| **AA ultra-low** | below ~125 MeV/n | BME | — | fixed in FLUKA |

chromo's FLUKA build does **not** link a UHE model (EPOS/SIBYLL/QGSJET). The class-level
ceiling `_max_sqrt_s_nn = 500 TeV` (in the center-of-mass frame) prevents requests above
the DPMJET-safe range.

---

## Installation

See [Building with FLUKA](../getting-started/building-with-fluka.md) for full instructions.
In brief:

```bash
export FLUPRO=$HOME/devel/FLUKA
bash scripts/install_fluka.sh
pip install --no-build-isolation -v -e .[test]
```

---

## InteractionType enum

The `InteractionType` enum controls which hadronic channels are active. The flag is a
three-digit decimal number where:

- **Units digit** — inelastic channel
- **Tens digit** — elastic channel
- **Hundreds digit** — electromagnetic dissociation (EMD)

| Value | Name | Channels | Notes |
|-------|------|----------|-------|
| `1` | `INELASTIC` | Inelastic only | Default |
| `10` | `ELASTIC` | Elastic only | |
| `11` | `INELA_ELA` | Inelastic + elastic | |
| `100` | `EMD` | EMD only | **Cross sections only** — see warning below |
| `101` | `INELA_EMD` | Inelastic + EMD | Recommended when EMD events are needed |
| `110` | `ELA_EMD` | Elastic + EMD | |
| `111` | `INELA_ELA_EMD` | All channels | |

!!! danger "`InteractionType.EMD` is for cross sections only"
    Requesting event generation with `interaction_type=InteractionType.EMD` (EMD-only,
    no inelastic channel) will abort FLUKA at runtime. For event generation that includes
    EMD, use `INELA_EMD` (101) or `INELA_ELA_EMD` (111) instead.

---

## Usage

```python
from chromo.models import Fluka
from chromo.models.fluka import InteractionType
from chromo.kinematics import FixedTarget

gen = Fluka(
    FixedTarget(100, "p", "O16"),
    interaction_type=InteractionType.INELA_EMD,
    seed=42,
)
for event in gen(10):
    print(event.final_state().pid)
```

### Cross sections

```python
xs = gen.cross_section()
print(f"σ_inel = {xs.inelastic:.2f} mb")
print(f"σ_emd  = {xs.emd:.2f} mb")
```

`cross_section()` works up to **1 PeV/nucleon** lab kinetic energy for both hadron and
photon projectiles. For nuclear targets, the inelastic cross section is also returned as
the production cross section (`xs.prod`).

### Constructor parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `evt_kin` | kinematics | — | Initial kinematics; target must be registered |
| `seed` | int or None | None | Random seed for FLUKA's Ranmar generator |
| `targets` | iterable | None | Extra nuclei to register beyond the 9 defaults |
| `interaction_type` | `InteractionType` | `INELASTIC` | Active hadronic channels |
| `transition_peanut_dpmjet` | float or None | 20 TeV | hA PEANUT↔DPMJET-3 threshold (GeV, lab kinetic) |
| `transition_peanut_dpmjet_smearing` | float or None | 0.0 | Smearing (±GeV) around the hA transition |
| `transition_rqmd_dpmjet` | float or None | 12.5 GeV/n | AA rQMD↔DPMJET-3 threshold (GeV/n, lab kinetic) |
| `transition_rqmd_dpmjet_smearing` | float or None | 2.0 GeV/n | Smearing (±GeV/n) around the AA transition |
| `enable_quasielastic` | bool | False | Enable quasi-elastic scattering |
| `rng_state_file` | Path or str | None (temp file) | File for persisting Ranmar RNG state across runs |

---

## Supported projectiles and targets

### Projectiles

Standard hadrons plus photons and light nuclei are supported for event generation.
Heavy-ion projectiles (A > 4) are accepted for cross-section queries only.

| Category | Examples |
|----------|---------|
| Hadrons | p, n, π±, K±, K⁰_S, K⁰_L |
| Photon | γ (`"gamma"`) |
| Light nuclei (cross section only for events) | He4, C12, ... |

### Targets

The following 9 nuclei are pre-registered by default at construction time. One additional
target slot is available (see [Hard material cap](#hard-material-cap-of-10-entries) below).

| Pre-registered targets |
|------------------------|
| p (free proton) |
| He4 |
| C12 |
| N14 |
| O16 |
| Ar40 |
| Fe56 |
| Cu63 |
| Pb208 |

To use a nucleus not in this list, pass it via the `targets` keyword:

```python
gen = Fluka(FixedTarget(100, "p", "Si28"), targets=["Si28"], seed=42)
```

---

## Limitations and caveats

### Single instantiation per process

Like all Fortran-based models in chromo, FLUKA uses global Fortran state and can only be
instantiated **once per Python process**. A second instantiation attempt will fail.

Tests work around this by running each model in a separate subprocess via
`tests/util.py::run_in_separate_process`.

---

### Hard material cap of 10 entries

!!! warning "Maximum 10 materials"
    FLUKA's `stpxyz.f` allocates its geometry with a compile-time `MEDFLK` upper bound
    of **10 regions**. This limit cannot be raised from Python or F2PY — it would require
    patching FLUKA's source.

The 9 default entries (`_DEFAULT_MATERIALS`) leave **1 free slot** for an extra target.
To use a nucleus outside the default list, replace elements by passing `targets=`:

```python
# Use Si28 as target — it occupies the 10th slot
gen = Fluka(FixedTarget(100, "p", "Si28"), targets=["Si28"], seed=42)
```

Exceeding 10 total unique materials raises `ValueError` at construction time. Adding new
targets at runtime is **not supported** — FLUKA's `STPXYZ` aborts on a second call, and
the underlying nuclear-table initialisation routines (`SETITB`, `DFATWG`) are not
accessible from Python.

---

### Energy ceilings differ between cross-section and event scopes

!!! warning "Different energy limits for `cross_section()` vs event generation"
    The usable energy range is **not the same** for cross-section queries and event
    generation:

| Scope | Projectile | Upper limit | Lower limit |
|-------|-----------|------------|-------------|
| `cross_section()` | hadron | ~1 PeV/nucleon lab | ~1 MeV/nucleon |
| `cross_section()` | photon | ~1 PeV/nucleon lab | ~1 MeV/nucleon |
| event generation | hadron (EVTXYZ) | **~20 TeV/nucleon lab** | ~1 MeV/nucleon |
| event generation | photon (PHNEVT) | ~100 TeV/nucleon lab | ~1 MeV/nucleon |

FLUKA's hadron event generator (EVTXYZ) crashes between 20 and 25 TeV/nucleon lab kinetic
energy. Stay below 20 TeV/nucleon for safe hadron event generation.

---

### Nuclear projectiles cannot generate events

!!! danger "AA event generation is not supported"
    `cross_section()` works for AA (nucleus–nucleus) kinematics, but calling `_generate()`
    with a heavy-ion projectile (A > 4) aborts inside FLUKA's EVTXYZ.

Use hadronic projectiles (p, π, K, n, γ) for event generation. If you need AA-like
final states, shoot hadrons at nuclear targets.

---

### EMD-only event generation aborts FLUKA

!!! danger "`EMD` interaction type is for cross sections only"
    `InteractionType.EMD` (value 100) activates only the electromagnetic dissociation
    channel. FLUKA will abort if you attempt event generation with this flag set and no
    inelastic channel enabled.

For event generation that includes EMD, use:

- `InteractionType.INELA_EMD` (101) — inelastic + EMD
- `InteractionType.INELA_ELA_EMD` (111) — all channels

---

### EMD cross section is zero for single-proton projectiles

EMD is driven by the Z²-enhanced Coulomb field of the projectile nucleus. For a bare
proton (Z = 1), the EMD cross section is physics-zero. Meaningful EMD cross sections
require nuclear projectiles (AA collisions).

---

### No beam records in HEPEVT

!!! note "HEPEVT contains only final-state particles"
    FLUKA's `FLLHEP` routine populates HEPEVT with **GENSTK ejectiles** and **RESNUC
    residual nuclei** only. There are no beam-particle entries at positions 0 and 1.

The generic `test_models_beam[Fluka]` test is marked `xfail` for this reason. Code that
assumes beam particles at the start of HEPEVT will not work with FLUKA events.

---

### `_set_stable` is a no-op

FLUKA's decay model is configured globally at initialisation time and cannot be changed
at runtime from Python. Calls to `_set_stable()` are silently ignored with an info-level
log message.

---

### e+/e- projectiles not supported

Electron and positron projectiles are not supported. Use `"gamma"` directly for
photon-induced interactions.

---

### RNG reproducibility requires Pythia8DecayHandler off

By default, chromo activates a `Pythia8DecayHandler` inside the FLUKA generator for
secondary decays. Pythia8 maintains its own internal RNG that is **independent from
FLUKA's Ranmar state** and is not seeded deterministically across processes. This means
that even with the same `seed=`, two runs may differ at the decay level.

For fully reproducible event records, disable the decay handler:

```python
gen = Fluka(FixedTarget(100, "p", "Pb208"), seed=42)
gen._activate_decay_handler(on=False)  # disable Pythia8 secondary decays
```

The Ranmar state can be saved and restored explicitly for checkpointing:

```python
gen.save_rng_state("checkpoint.dat")
# ... later ...
gen.load_rng_state("checkpoint.dat")
```

---

## CompositeTarget support

FLUKA supports `CompositeTarget` (e.g. Air = N + O + Ar). All components must be
pre-registered materials. The default material list already includes N14, O16, and Ar40,
so Air works out of the box:

```python
from chromo.util import CompositeTarget
from chromo.kinematics import FixedTarget

air = CompositeTarget([("N14", 0.78), ("O16", 0.21), ("Ar40", 0.01)])
gen = Fluka(FixedTarget(1000, "p", air), seed=42)
```

---

## See also

- [Building with FLUKA](../getting-started/building-with-fluka.md) — installation instructions
- [Model overview](overview.md) — comparison with other chromo backends
