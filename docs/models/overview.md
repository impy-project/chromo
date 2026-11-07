# Model Overview

chromo wraps 9 model families with 26 model classes. This page summarizes the capabilities and limitations of each.

!!! note "Citing models"
    chromo only provides a Python interface. The physics models themselves are developed independently. When using a model in public work, always cite the original authors. See [Citation](../citation.md).

## Capability Table

| Model Class | Family | Projectiles | Targets | E_cm min | Platform | Notes |
|-------------|--------|-------------|---------|----------|----------|-------|
| `Sibyll21` | SIBYLL | h | N, A (A<=20) | 10 GeV | all | Legacy |
| `Sibyll23c` | SIBYLL | h | N, A (A<=20) | 10 GeV | all | |
| `Sibyll23d` | SIBYLL | h | N, A (A<=20) | 10 GeV | all | Recommended |
| `Sibyll23e` | SIBYLL | h | N, A (A<=20) | 10 GeV | all | Latest tune |
| `Sibyll23dStarMixed` | SIBYLL | h, A (A<=56) | N, A (A<=20) | 10 GeV | all | Nuclear proj |
| `Sibyll23eStarMixed` | SIBYLL | h, A (A<=56) | N, A (A<=20) | 10 GeV | all | Nuclear proj |
| `QGSJet01d` | QGSJet | h, A | N, A | 10 GeV | all | Legacy |
| `QGSJetII03` | QGSJet | h, A | N, A | 10 GeV | all | |
| `QGSJetII04` | QGSJet | h, A | N, A | 10 GeV | all | Widely used |
| `QGSJetIII` | QGSJet | h, A | N, A | 10 GeV | all | Latest |
| `EposLHC` | EPOS | h, A | N, A | 6 GeV | all | |
| `EposLHCR` | EPOS | h, A | N, A | 6 GeV | all | Updated tune |
| `EposLHCRHadrRescattering` | EPOS | h, A | N, A | 6 GeV | all | + UrQMD rescattering (slow) |
| `DpmjetIII307` | DPMJET | h, A (A<=280) | A | 1 GeV | all | Legacy |
| `DpmjetIII191` | DPMJET | h, A (A<=280) | A | 1 GeV | all | Modern |
| `DpmjetIII193` | DPMJET | h, A (A<=280) | A | 1 GeV | all | Latest |
| `Phojet112` | PHOJET | h | N | 10 GeV | all | Legacy |
| `Phojet191` | PHOJET | h, gamma | N, gamma | 10 GeV | all | |
| `Phojet193` | PHOJET | h, gamma | N, gamma | 10 GeV | all | |
| `Pythia6` | PYTHIA | h | N | 10 GeV | all | |
| `Pythia8` | PYTHIA 8 | h, e, gamma | N, e, gamma | 10 GeV | Linux, macOS | See [guide](pythia8-guide.md) |
| `Pythia8Cascade` | PYTHIA 8 | h, A | A (A>1) | 10 GeV | Linux, macOS | h+A single coll. |
| `Pythia8Angantyr` | PYTHIA 8 | h, A | A | 20 GeV | Linux, macOS | Heavy-ion Glauber |
| `Sophia20` | SOPHIA | gamma | N | -- | all | Photoproduction only |
| `UrQMD34` | UrQMD | h, A | A | 2 GeV | Linux, macOS | Transport model |
| `Fluka` | FLUKA | h, gamma, A | A | ~1 MeV/n | Linux, macOS | License required |

**Legend:** *h* = hadrons (p, n, pi+/-, K+/-, ...), *N* = nucleon (p or n), *A* = nucleus, *gamma* = photon, *e* = electron/positron

## Choosing a Model

**Cosmic-ray air showers:** SIBYLL 2.3d/e, EPOS-LHC-R, or QGSJet-II-04/III are the standard choices. All handle h+Air via `CompositeTarget`.

**Heavy-ion collisions:** DPMJET-III 19.x, EPOS-LHC, Pythia8 Angantyr, or UrQMD. For Pythia 8, see the [Pythia 8 Guide](pythia8-guide.md).

**Photon interactions:** SOPHIA (gamma+N), PHOJET 19.x (gamma+gamma, gamma+N), or Pythia 6/8 (gamma+N, gamma+gamma).

**Low-energy nuclear:** FLUKA (down to ~1 MeV/nucleon) or UrQMD (from 2 GeV CMS).

**Collider physics (pp):** Pythia 8 for the most detailed simulation, or any of the cosmic-ray models for comparison.
