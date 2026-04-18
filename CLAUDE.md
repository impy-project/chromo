# CLAUDE.md

## Project Overview

**chromo** (Cosmic Ray and HadROnic interactiOn MOnte-carlo frontend) is a Python package providing a unified interface to particle physics event generators (SIBYLL, QGSJet, PYTHIA, DPMJET, EPOS, UrQMD, SOPHIA). It wraps Fortran/C/C++ backends via F2PY and pybind11.

## Build System

The build uses **Meson** (via `mesonpy`). Key files: `meson.build`, `pyproject.toml`.

```bash
# Install from source (editable, with test deps)
pip install meson-python numpy ninja pre-commit
pip install --no-build-isolation -v -e .[test]

# Requires: C/C++ compiler, Fortran compiler (gcc/gfortran), Python 3.9+
# Submodules must be checked out: git submodule update --init --recursive
```

Compiled modules go to `build/cp*/`. Model enable/disable is controlled in `pyproject.toml` under `[tool.chromo]`.

## Running Tests

```bash
# Download required data files before running tests
python scripts/download_data.py

# Run all tests (parallel by default via pytest-xdist)
# This dev box has 32-40 cores available; use -n 32. CI runners use 2-3.
python -m pytest -vv -n 32 # in CI/cloud use 2 or 3 processes.

# Run tests serially
python -m pytest -n 0

# Skip slow generator tests
python -m pytest -k "not test_generators"

# Debug a model in current process (only works once per model)
DEBUG=10 python -m pytest -n 0 -k test_name
```

Data files are cached in `~/.cache/chromo`. Download is not thread-safe, so run `scripts/download_data.py` before parallel test execution.

## Code Style

Styling, linting and formatting is configured via pre-commit hooks and `pyproject.toml`:
Run `pre-commit run -a` to apply all hooks to the codebase, and fix errors before committing.

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
Each model file (e.g., `sibyll.py`, `qgsjet.py`, `pythia8.py`) subclasses `MCRun`/`MCEvent` and interfaces with compiled backends. Models can only be instantiated **once per process** due to Fortran global state — tests run each model in a separate subprocess. The relevant model classes are exposed in `src/chromo/models/__init__.py:__all__` for user access, extra models can be accessed through `src/chromo/models/_extra_models.py` if needed, whereas those with "_DEV" suffix are for development and testing only and are not located in the source code base.

### Pythia8 Nuclear Modes (`src/chromo/models/pythia8.py`)

Three Pythia8 model classes with distinct physics scopes:

- **`Pythia8`** — Standard Pythia8 for hN, ee, γγ, γN collisions. Extended projectiles include strange/charm/bottom hadrons. No nuclear targets.
- **`Pythia8Cascade`** — PythiaCascade plugin for single-collision h+A mode. Nuclear projectiles decomposed into Z protons + (A-Z) neutrons. Targets: nuclei with A>1. Uses `slowDecays=True` (cosmic-ray convention). RNG state save/restore via both internal Pythia instances (pythiaMain + pythiaColl). Supports `CompositeTarget` via `_composite_plan`.
- **`Pythia8Angantyr`** — Glauber heavy-ion model for hA/AA with precomputed tables (20 GeV–20 PeV CMS, `_ecm_min=20 GeV`). Targets: nuclei only (no proton/neutron). Live target switching via `setBeamIDs`. Cross sections: `cross_section()` returns fast parametric estimates (PythiaCascade `nCollAvg` formula); `glauber_cross_section(n_trials)` runs GlauberOnly MC for precise Angantyr values. Supports `CompositeTarget`.

C++ bindings in `src/cpp/_pythia8.cpp`: `PythiaCascadeForChromo` wraps `PythiaCascade` with numpy array extraction, `set_may_decay` on both internal Pythia instances, and `getRndmState`/`setRndmState` for RNG state serialization.

### Fortran/C++ Sources
- `src/fortran/` — Fortran source for SIBYLL, QGSJet, DPMJET, EPOS, UrQMD, SOPHIA, PYTHIA6
- `src/cpp/` — PYTHIA 8 C++ source with pybind11 bindings
- `scripts/generate_f2py.py` — Custom F2PY wrapper generation

### I/O (`src/chromo/writer.py`)
Output formats: HepMC3 (text/gzip), ROOT (via uproot), SVG visualization.

### CLI (`src/chromo/cli.py`)
Command-line interface compatible with CRMC. Entry point: `chromo` command.

## Windows Build Support

Windows builds are not yet **experimentally supported**. Key details:

### Windows-Specific Build Requirements
- **Compiler:** MinGW-w64 with gfortran (GCC 13.x recommended)
- **Command line length limitation:** Windows cmd.exe has an 8,191 character limit. The build system uses response files (`--source-file-list`) to pass source file lists to `generate_f2py.py`, avoiding this limit.

### Response File Implementation
- `scripts/generate_f2py.py` supports `--source-file-list <file>` to read source paths from a text file
- `meson.build` creates `.txt` files with source lists for all models (especially critical for DPMJET 19.1/19.3 with 700+ files)
- This reduces command lines from 42,000+ characters to ~200 characters

### Known Windows Limitations
- **PYTHIA8:** Not built on Windows, yet. Conditionally skipped in `meson.build` line 356: `if 'pythia8' in enabled_models and host_machine.system() != 'windows'`
- **UrQMD:** May fail to compile with gfortran due to non-standard Fortran extensions. If so, disable in `pyproject.toml`.
- **QGSJet set_stable():** Known issue on Windows, tests marked xfail in `tests/test_setstable.py`.
- **Test parallelism:** Use `-n 2` instead of `-n 3` on Windows for stability.

### Windows Testing Command
```bash
# Skip PYTHIA8 and problematic tests
python -m pytest -vv -k "not (Pythia8 or test_decay_handler)" -n 2
```

## FLUKA integration

FLUKA 2025.1 is an optional, license-restricted backend. Not built in
public CI.

### Install

```bash
export FLUPRO=$HOME/devel/FLUKA
export FLUKA_ARCHIVE_DIR=$HOME/devel/FLUKA-dev  # directory with the .tar.gz archives
bash scripts/install_fluka.sh
```

Persist `export FLUPRO=$HOME/devel/FLUKA` in your shell rc.

### Build chromo with Fluka

```bash
pip install --no-build-isolation -v -e .[test]
```

The meson `fluka` block fails fast if `$FLUPRO` is unset or any of the
required archives (`libflukahp.a`, `libdpmmvax.a`, `librqmdmvax.a`,
`latestRQMD/librqmd.a`, `interface/libdpmjet*.a`, `interface/dpmvers`)
are missing.

### Usage

```python
from chromo.models import Fluka
from chromo.models.fluka import InteractionType
from chromo.kinematics import FixedTarget

gen = Fluka(FixedTarget(100, "p", "O16"),
            interaction_type=InteractionType.INELA_EMD,
            seed=42)
for event in gen(10):
    print(event.final_state().pid)
```

### Caveats & FLUKA-internal limitations

- **Single instantiation per Python process.** Fortran globals; tests use
  `tests/util.py::run_in_separate_process`.
- **Hard material cap of 10 entries.** FLUKA's shipped `stpxyz.f:256` has
  a compile-time `MEDFLK` upper bound of 10 (empirically confirmed: index
  11 fires a Fortran runtime array-bound error). Cannot be raised from
  Python/F2PY; would require patching FLUKA's source. Chromo's default
  `_DEFAULT_MATERIALS` list has 9 entries (`p`, `He4`, `C12`, `N14`,
  `O16`, `Ar40`, `Fe56`, `Cu63`, `Pb208`) leaving 1 slot for an extra
  target. Use `targets=["Si28", ...]` to swap elements; exceeding 10
  raises `ValueError`.
- **Energy ceiling at 300 TeV CMS** (`_ecm_max`). `PPTMAX` is set to
  the construction-kinematics `plab`, so DPMJET Glauber tables cover
  the full requested range. `cross_section()` works at all energies
  below the ceiling.
- **Light-nucleus projectiles (d, t, 3He, 4He) work** for both cross
  sections and event generation via dedicated FLUKA PAPROP codes.
- **Heavy-ion projectiles (A > 4) are not yet supported** for event
  generation (pending upstream support). `cross_section()` works for
  AA kinematics via SGMXYZ.
- **EMD-only event generation aborts FLUKA.** Use `InteractionType.EMD`
  for cross-section queries only. For event generation, combine with
  `INELA_EMD` (101) or `INELA_ELA_EMD` (111).
- **EMD cross section is zero for single-proton projectiles** (physics:
  needs a Z²-enhanced field, so meaningful EMD only appears for AA).
- **No beam records in HEPEVT.** FLUKA's `FLLHEP` populates HEPEVT with
  GENSTK ejectiles and RESNUC residuals. `FlukaEvent._prepend_initial_beam`
  is intentionally a no-op. The generic `test_models_beam[Fluka]` is
  xfailed for this reason (see `tests/test_common.py`).
- **`_set_stable` is a no-op.** FLUKA's decay model is global and not
  runtime-configurable.
- **`e+/e-` projectiles are not yet supported** (SGMXYZ returns xsec
  but EVTXYZ aborts on some targets). Pending upstream clarification.
- **RNG reproducibility requires the Pythia8DecayHandler off.** Pythia8's
  own RNG is independent from FLUKA's Ranmar state and isn't seeded
  deterministically across processes. For fully reproducible event
  records: `Fluka(... )._activate_decay_handler(on=False)`.

### Disable

Remove `"fluka"` from `[tool.chromo] enabled-models` in `pyproject.toml`.

## Data Files (iamdata)

Each model downloads a versioned zip from GitHub releases into `src/chromo/iamdata/<model>/`.
The mechanism lives in `src/chromo/util.py:_cached_data_dir(url)`:

1. Derives `model_name` and `vname_stem` (e.g. `Pythia8_v006`) from the URL.
2. If `iamdata/<model>/<vname_stem>` (a version marker file) already exists → returns immediately, no download.
3. Otherwise clears the model dir, copies the zip from `~/.cache/chromo/` (CI) or downloads it, extracts it, then creates the version marker file.

**To bump a model's data version** (e.g. Pythia8 v005 → v006):
1. Update `_version` and `_data_url` in the model file (`src/chromo/models/<model>.py`).
2. Repopulate `iamdata/<model>/` from the submodule (xmldoc, tunes, pdfdata for Pythia8).
3. Remove the old version marker (`Pythia8_v005`) and create the new one (`Pythia8_v006`) — just a plain file, content doesn't matter.
4. Create the zip: `cd src/chromo/iamdata && zip -r <Model>_v00N.zip <Model>/xmldoc <Model>/pdfdata <Model>/tunes`
   The zip must contain `<Model>/subdir/...` paths (NOT the version file — that is created after extraction).
5. Upload zip to `https://github.com/impy-project/chromo/releases/tag/zipped_data_v1.0`.
6. Add the new version name to the cache download loop in `.github/workflows/test.yml` and `build.yml`, and also bump/rotate the GitHub Actions cache key there whenever the cached zip list changes (cache entries are immutable for a given key).

Note: `iamdata/` is git-ignored; only the model `.py` and CI workflow files are committed.

## Important Constraints

- **Single instantiation:** Most Fortran models use global state and can only be initialized once per Python process. Tests work around this by running models in separate subprocesses.
- **Symbol management:** Exported symbols are defined per model in `meson.build`. Tests in `test_exposed_symbols.py` verify these. Keep `meson.build` function lists and test expectations in sync.
- **Adding a new model:** Update `meson.build` with sources and symbol lists, add to `[tool.chromo] enabled-models` in `pyproject.toml`, create model wrapper in `src/chromo/models/`, add tests.
- **Windows command line limits:** When adding models with many source files, ensure they use the response file pattern in `meson.build` (see DPMJET section for example).
- **Debugging Python code inline:** When you need to run short Python snippets, write them to a temp.py file and execute that instead of using `python -c`. This avoids permission issues.
