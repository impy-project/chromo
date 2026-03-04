# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**chromo** (Cosmic Ray and HadROnic interactiOn MOnte-carlo frontend) is a Python package providing a unified interface to particle physics event generators (SIBYLL, QGSJet, PYTHIA, DPMJET, EPOS, UrQMD, SOPHIA). It wraps Fortran/C/C++ backends via F2PY and pybind11.

## Build System

The build uses **Meson** (via `mesonpy`). Key files: `meson.build`, `pyproject.toml`.

```bash
# Install from source (editable, with test deps)
pip install --prefer-binary -v -e .[test]

# On Windows, set these env vars:
#   CMAKE_GENERATOR="MinGW Makefiles"
#   FC=gfortran

# Requires: C/C++ compiler, Fortran compiler (gcc/gfortran), Python 3.9+
# Submodules must be checked out: git submodule update --init --recursive
```

Compiled modules go to `build/cp*/`. Model enable/disable is controlled in `pyproject.toml` under `[tool.chromo]`.

## Running Tests

```bash
# Run all tests (parallel by default via pytest-xdist)
python -m pytest -vv -n 16

# Run tests serially
python -m pytest -n 0

# Skip slow generator tests
python -m pytest -k "not test_generators"

# Debug a model in current process (only works once per model)
DEBUG=10 python -m pytest -n 0 -k test_name

# Download required data files before running tests
python scripts/download_data.py
```

Data files are cached in `~/.cache/chromo`. Download is not thread-safe, so run `scripts/download_data.py` before parallel test execution.

## Code Style

- **Formatter:** ruff format (line length 88)
- **Linter:** ruff (configured in `pyproject.toml`)
- **Docstrings:** NumPy convention
- **Pre-commit hooks:** configured in `.pre-commit-config.yaml` (ruff, trailing whitespace, etc.)

## Architecture

### Core Classes (`src/chromo/common.py`)
- **`MCRun`** — Abstract base for all generators. Subclasses implement `_generate()`, `_set_kinematics()`, `_cross_section()`, `_set_stable()`.
- **`MCEvent`** — Base event class wrapping the HEPEVT common block from Fortran.
- **`EventData`** — Pickleable dataclass holding particle arrays (pid, px, py, pz, en, mass, vertices, family relations). Provides filtering (`final_state()`, `final_state_charged()`), frame transformations, and HepMC export.
- **`CrossSectionData`** — Cross-section results (total, elastic, inelastic, diffractive).

### Kinematics (`src/chromo/kinematics.py`)
- **`CenterOfMass`** / **`FixedTarget`** — Specify collision frame. Handles energy/momentum conversions and CMS ↔ lab frame transformations.
- **`CompositeTarget`** — Mixed materials (e.g., Air = N + O + Ar).

### Model Implementations (`src/chromo/models/`)
Each model file (e.g., `sibyll.py`, `qgsjet.py`, `pythia8.py`) subclasses `MCRun`/`MCEvent` and interfaces with compiled backends. Models can only be instantiated **once per process** due to Fortran global state — tests run each model in a separate subprocess.

### Fortran/C++ Sources
- `src/fortran/` — Fortran source for SIBYLL, QGSJet, DPMJET, EPOS, UrQMD, SOPHIA, PYTHIA6
- `src/cpp/` — PYTHIA 8 C++ source with pybind11 bindings
- `scripts/generate_f2py.py` — Custom F2PY wrapper generation

### I/O (`src/chromo/writer.py`)
Output formats: HepMC3 (text/gzip), ROOT (via uproot), SVG visualization.

### CLI (`src/chromo/cli.py`)
Command-line interface compatible with CRMC. Entry point: `chromo` command.

## Important Constraints

- **Single instantiation:** Most Fortran models use global state and can only be initialized once per Python process. Tests work around this by running models in separate subprocesses.
- **Symbol management:** Exported symbols are defined per model in `meson.build`. Tests in `test_exposed_symbols.py` verify these. Keep `meson.build` function lists and test expectations in sync.
- **Adding a new model:** Update `meson.build` with sources and symbol lists, add to `[tool.chromo] enabled-models` in `pyproject.toml`, create model wrapper in `src/chromo/models/`, add tests.
