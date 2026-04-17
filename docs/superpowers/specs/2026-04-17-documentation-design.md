# chromo Documentation Design Spec

**Date:** 2026-04-17
**Status:** Approved
**Branch:** `docs/mkdocs` (off `feature_fluka`)

## Overview

Create a user-facing documentation site for chromo using MkDocs + Material theme. Primary audience: physicists who `pip install chromo` and want to simulate particle interactions. The docs should explain usage, expose model limitations, embed existing Jupyter notebooks, and auto-generate API reference from docstrings.

## Tooling

- **Static site generator:** MkDocs with Material for MkDocs theme
- **API docs:** `mkdocstrings[python]` (auto-generated from docstrings)
- **Notebook rendering:** `mkdocs-jupyter` (embed existing notebooks as doc pages)
- **Deployment:** GitHub Pages (via `mkdocs gh-deploy` or CI action)
- **Python deps:** `mkdocs-material`, `mkdocstrings[python]`, `mkdocs-jupyter`

## Navigation Structure

Task-first user guide with an architecture reference page:

```
chromo Documentation
├── Home (index)
├── Getting Started
│   ├── Installation
│   ├── Your First Simulation
│   └── Building with FLUKA
├── User Guide
│   ├── Kinematics
│   ├── Running a Generator
│   ├── Working with Events
│   ├── Cross Sections
│   ├── Output Formats
│   └── Decay Handling
├── Models
│   ├── Model Overview (capability table)
│   ├── Pythia 8 Guide (decision tree: Pythia8 vs Cascade vs Angantyr)
│   ├── SIBYLL Variants (2.1/2.3c/d/e/Star)
│   ├── DPMJET & PHOJET (3.0.7 vs 19.x, relationship)
│   └── FLUKA (caveats, InteractionType, limitations)
├── Architecture (paper diagram, layered design)
├── Examples (embedded notebooks)
├── API Reference (auto-generated)
│   ├── chromo.kinematics
│   ├── chromo.common
│   ├── chromo.models
│   ├── chromo.writer
│   └── chromo.cli
└── Citation
```

## File Layout

```
docs/
├── index.md
├── getting-started/
│   ├── installation.md
│   ├── first-simulation.md
│   └── building-with-fluka.md
├── user-guide/
│   ├── kinematics.md
│   ├── running-a-generator.md
│   ├── working-with-events.md
│   ├── cross-sections.md
│   ├── output-formats.md
│   └── decay-handling.md
├── models/
│   ├── overview.md
│   ├── pythia8-guide.md
│   ├── sibyll-variants.md
│   ├── dpmjet-phojet.md
│   └── fluka.md
├── architecture.md
├── examples/                # references to ../examples/*.ipynb
├── api/
│   ├── kinematics.md
│   ├── common.md
│   ├── models.md
│   ├── writer.md
│   └── cli.md
└── citation.md
mkdocs.yml
```

## Page Content Details

### Home (index.md)

Project tagline from README, the chromo SVG logo, a 5-line "hello world" example, links to installation and user guide. Short and inviting.

### Getting Started

**Installation:** Three paths: `pip install chromo` (recommended), from-source with meson (existing `doc/dev_docs.md` content), Docker. Platform notes (macOS Xcode issue, Windows limitations).

**Your First Simulation:** Walk through the README example with explanations of each step: what `CenterOfMass` does, what the generator iterator returns, what `final_state_charged()` filters. Links to "where to go next."

**Building with FLUKA:** FLUKA is license-restricted. Prerequisites: obtain FLUKA archives, run `install_fluka.sh`, set `$FLUPRO`, build chromo from source. Clear checklist format.

### User Guide

**Kinematics:** `CenterOfMass`, `FixedTarget`, `CompositeTarget` explained with examples. Energy units (`GeV`, `TeV`, `PeV`, `EeV`). How frame conversions work.

**Running a Generator:** The critical **single-instantiation constraint** due to Fortran global state. Why you can only create one model per process. How to use `multiprocessing` or subprocesses to work around it. The `generator(n_events)` iteration protocol. Code examples.

**Working with Events:** `EventData` fields (pid, px, py, pz, en, mass, status, vertices). Filtering methods (`final_state()`, `final_state_charged()`). Frame transforms. Practical "how do I get pT of charged pions" patterns.

**Cross Sections:** `CrossSectionData` fields and what each means (total, inelastic, elastic, diffractive components, EMD). Which models fill which fields.

**Output Formats:** HepMC3 (text/gzip), ROOT (via uproot), SVG visualization. CLI usage for output generation.

**Decay Handling:** `Pythia8DecayHandler`, how to control which particles decay, interaction with model-native decay settings.

### Models

**Model Overview:** Capability table covering all 26 model classes. Columns: Model, Versions, Projectiles, Targets, Energy Range, Platform, Notes. Derived from `_projectiles`, `_targets`, `_ecm_min` class attributes but hand-written for clarity.

**Pythia 8 Guide:** Decision tree:
- `Pythia8` -- standard hN, ee, gamma. No nuclear targets. Extended projectile set (strange/charm/bottom hadrons).
- `Pythia8Cascade` -- PythiaCascade plugin for single-collision h+A. Nuclear projectiles decomposed into nucleons. Targets: A>1. `slowDecays=True`. Supports `CompositeTarget`.
- `Pythia8Angantyr` -- Glauber heavy-ion model for hA/AA. Precomputed tables 20 GeV-20 PeV CMS. No proton/neutron targets. `cross_section()` vs `glauber_cross_section()`. Supports `CompositeTarget`.

Each with code example, supported kinematics, and gotchas.

**SIBYLL Variants:** Table of 2.1 / 2.3c / 2.3d / 2.3e / Star. What changed between versions. Target limits (A<=20 classic, A<=56 Star). Guidance on which to use for which purpose.

**DPMJET & PHOJET:** Relationship (PHOJET = photon interaction mode of DPMJET). Version matrix (3.0.7 legacy, 19.1/19.3 modern). Nuclear capabilities (A up to 280). Extended projectiles. Recommendations.

**FLUKA:** InteractionType enum, 10-material cap, energy ceilings (cross-section up to 1 PeV/nucleon vs event generation up to ~20 TeV/nucleon for hadrons), nuclear projectile limitation (cross sections yes, events no), EMD caveats, seed reproducibility, no beam records in HEPEVT, `_set_stable` is a no-op. User-friendly version of CLAUDE.md caveats.

### Architecture

Reproduce the paper's program structure diagram (embedded image or Mermaid). Explain the three layers:
1. `chromo.kinematics` -- user interface unification (frame conversions, particle specification)
2. `chromo.common` -- wrapping and technical IO unification (MCRun, MCEvent, EventData)
3. `chromo.models.*` -- specific implementations for each generator family

Links to API reference for each component.

### Examples

Embedded notebooks via `mkdocs-jupyter`, ordered by topic:
1. `compare_models.ipynb` -- comparing generators
2. `cross_section.ipynb` -- querying cross sections
3. `nuclear_cross_sections.ipynb` -- nuclear/heavy-ion
4. `count_final_state_particles.ipynb` -- working with arrays
5. `gamma_p.ipynb` / `gamma_gamma_example.ipynb` -- photon interactions
6. `hepmc_io_and_visualization.ipynb` -- output and visualization
7. `hyperon_feed_down.ipynb` -- advanced physics
8. `decayhandler.ipynb` -- controlling decays
9. `write_root.ipynb` -- ROOT output
10. `fluka_dpmjet_residual_nuclei.ipynb` -- FLUKA-specific

Each with a short intro sentence.

### API Reference

Auto-generated with `mkdocstrings-python`:
- `chromo.kinematics` -- `CenterOfMass`, `FixedTarget`, `CompositeTarget`, `EventKinematicsBase`
- `chromo.common` -- `MCRun`, `MCEvent`, `EventData`, `CrossSectionData`
- `chromo.models` -- one entry per model class from `__all__`
- `chromo.writer` -- output writers
- `chromo.cli` -- CLI reference

### Citation

BibTeX from README. Table linking each model to its INSPIRE record.

## Phased Implementation

**Phase 1:** Scaffold `mkdocs.yml`, write user guide pages (getting started, kinematics, events, model overview table, output), embed key notebooks, set up `mkdocstrings` for core classes.

**Phase 2:** Add deep-dive model pages (Pythia8 guide, SIBYLL variants, DPMJET & PHOJET), architecture page with paper diagram.

**Phase 3:** FLUKA section (install, build-from-source, usage, caveats, example notebook).

## Branch Strategy

- Create `docs/mkdocs` branch off `feature_fluka`
- All docs work happens on `docs/mkdocs`
- When `feature_fluka` merges to `main`, rebase `docs/mkdocs` onto `main`
