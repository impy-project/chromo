# SIBYLL Variants

SIBYLL is a hadronic interaction model designed for cosmic-ray air shower simulations. chromo
ships six SIBYLL variants spanning two generations and two projectile-range families.

---

## Variant comparison table

| Class | Based on | Projectiles | Max target A | Notes |
|-------|----------|-------------|--------------|-------|
| `Sibyll21` | 2.1 | Standard hadrons | 20 | Legacy. No Ω⁻ / anti-Ω production. |
| `Sibyll23c` | 2.3c | Standard hadrons | 20 | Updated charm production tune. |
| `Sibyll23d` | 2.3d | Standard hadrons | 20 | **Recommended default** for air showers. |
| `Sibyll23e` | 2.3e | Standard hadrons | 20 | Latest tune; improved forward production. |
| `Sibyll23dStarMixed` | 2.3d + SIBYLL* | Hadrons + nuclei up to A=56 | 20 | Nuclear projectiles; same target limit. |
| `Sibyll23eStarMixed` | 2.3e + SIBYLL* | Hadrons + nuclei up to A=56 | 20 | Nuclear projectiles; same target limit. |

---

## When to use which variant

**`Sibyll23d`** — the standard choice for cosmic-ray air shower simulations. This version is used
in established benchmark studies and is well-validated against accelerator and air shower data.
Use it unless you have a specific reason to deviate.

**`Sibyll23e`** — the most recent tune, with improved modelling of forward particle production.
Choose this if you want the latest physics or are comparing to very recent data sets.

**Star variants (`Sibyll23dStarMixed`, `Sibyll23eStarMixed`)** — SIBYLL* extends the projectile
range to nuclear primaries up to A = 56. Use these only when you need to simulate nuclear
projectiles such as iron-56 or carbon-12 striking air nuclei. The physics of the hadronic
interaction is identical to the parent 2.3d/2.3e version; SIBYLL* adds the nuclear-projectile
decomposition on top.

**`Sibyll21`** — kept for backward compatibility and reproducibility of older results. Note that
it does not produce Ω⁻ and anti-Ω baryons, which matters if you need those species.

**`Sibyll23c`** — intermediate version with an updated charm tune. Useful for charm-physics
comparisons but superseded by 2.3d/2.3e for general use.

---

## Target limits

Classic SIBYLL variants (2.1, 2.3c/d/e) support targets up to **A = 20**, which covers the
main atmospheric nuclei: nitrogen-14 (N14), oxygen-16 (O16), and argon-40 (Ar40). Heavier
targets are not supported.

SIBYLL* variants keep the same target limit (A ≤ 20) but extend the **projectile** range to
A ≤ 56, enabling iron-56 (Fe56) and lighter nuclei as primaries.

!!! warning "Target mass limit"
    Requesting a target with A > 20 in any SIBYLL variant will be rejected by the kinematics
    validation layer. Use DPMJET or Pythia8Angantyr for heavier targets such as lead.

---

## Code examples

### Standard air shower simulation with Sibyll23d

```python
from chromo.models import Sibyll23d
from chromo.kinematics import FixedTarget

# 100 TeV proton on nitrogen-14 (typical cosmic-ray air shower target)
gen = Sibyll23d(FixedTarget(1e5, "p", "N14"), seed=1)

for event in gen(200):
    fs = event.final_state()
    print(f"N_particles = {len(fs)}, total E = {fs.en.sum():.1f} GeV")
```

### Nuclear projectile with a Star variant

```python
from chromo.models import Sibyll23dStarMixed
from chromo.kinematics import FixedTarget

# Iron-56 on oxygen-16 at 10 PeV/nucleon (cosmic-ray heavy primary)
gen = Sibyll23dStarMixed(FixedTarget(1e7, (26, 56), "O16"), seed=42)

for event in gen(50):
    fs = event.final_state_charged()
    print(f"Nch = {len(fs)}")
```

### Using a CompositeTarget for air

```python
from chromo.models import Sibyll23d
from chromo.kinematics import FixedTarget, CompositeTarget

air = CompositeTarget([("N14", 0.78), ("O16", 0.21), ("Ar40", 0.01)])
gen = Sibyll23d(FixedTarget(1e4, "p", air), seed=5)

for event in gen(100):
    print(event.final_state().pid[:10])
```

### Cross sections

```python
from chromo.models import Sibyll23d
from chromo.kinematics import FixedTarget

gen = Sibyll23d(FixedTarget(1e5, "p", "N14"), seed=1)
xs = gen.cross_section()
print(f"σ_total   = {xs.total:.2f} mb")
print(f"σ_inel    = {xs.inelastic:.2f} mb")
print(f"σ_elastic = {xs.elastic:.2f} mb")
```

---

## Physics notes

- All SIBYLL variants use a quark-gluon string (minijet + soft string) model for the hadronic
  interaction. The main differences between versions are tuning of parameters and production
  of secondary hadrons.
- Forward particle production (the "fragmentation" region) is especially important for air shower
  development. SIBYLL 2.3e contains the most refined treatment of this region.
- Charm production in SIBYLL 2.3c and later is relevant for the prompt atmospheric neutrino flux.
  Use `Sibyll23c` or later if charm is your observable of interest.
