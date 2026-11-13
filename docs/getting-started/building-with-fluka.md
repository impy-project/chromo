# Building with FLUKA

[FLUKA](https://fluka.cern) is an optional, license-restricted backend. It is not included in the standard chromo wheel and must be built from source.

## Prerequisites

1. **Obtain FLUKA 2025.1** -- FLUKA requires a license from [fluka.cern](https://fluka.cern). Download the two archives (the main FLUKA package and the data files).

2. **Place the archives** in a known directory (e.g., `$HOME/devel/FLUKA-dev/`).

3. **C/C++ and Fortran compilers** -- same as for building chromo from source.

## Install FLUKA

```bash
# Set the FLUKA installation directory
export FLUPRO=$HOME/devel/FLUKA

# Run the install script (expects archives in $HOME/devel/FLUKA-dev/)
# Set FLUKA_ARCHIVE_DIR to override the archive location
bash scripts/install_fluka.sh
```

Add `export FLUPRO=$HOME/devel/FLUKA` to your shell profile (`.bashrc`, `.zshrc`, etc.) so it persists across sessions.

## Build chromo with FLUKA

```bash
# From the chromo source directory
pip install --no-build-isolation -v -e .[test]
```

The Meson build system detects `$FLUPRO` and links against the FLUKA libraries. It will fail fast if `$FLUPRO` is unset or required archives are missing.

## Verify the Installation

```python
from chromo.models import Fluka
from chromo.kinematics import FixedTarget

gen = Fluka(FixedTarget(100, "p", "O16"), seed=42)
for event in gen(5):
    print(f"Event {event.nevent}: {len(event.final_state())} final-state particles")
```

## Disabling FLUKA

Remove `"fluka"` from the `[tool.chromo] enabled-models` list in `pyproject.toml`.
