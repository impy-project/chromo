# Installation

## From PyPI (recommended)

The easiest way to install chromo is from the pre-compiled binary wheel:

```bash
pip install chromo
```

Binary wheels are available for:

- **Python:** 3.9+
- **Platforms:** Linux (x86_64, aarch64), macOS (x86_64, arm64), Windows (x86_64)

## From Source

Building from source requires C/C++ and Fortran compilers (GCC + gfortran recommended).

```bash
# Clone with submodules
git clone --recursive https://github.com/impy-project/chromo
cd chromo

# Install build dependencies
pip install meson-python numpy ninja

# Install in editable mode
pip install --no-build-isolation -v -e .[test]
```

Check `[build-system.requires]` in `pyproject.toml` for the full list of build dependencies.

### macOS Notes

- Install GCC and gfortran via Homebrew: `brew install gcc`
- Xcode Command Line Tools version 14 has a known linker bug with C++ compiled by GCC. Workaround: downgrade to version 13.4 and disable automatic updates.

### Windows Notes

- Requires MinGW-w64 with gfortran (GCC 13.x recommended).
- PYTHIA 8 is not yet available on Windows.
- UrQMD may fail to compile due to non-standard Fortran extensions.
- Use `-n 2` (not `-n 3`) for parallel tests.

## From Source in Docker

For a verified build environment:

```bash
git clone --recursive https://github.com/impy-project/chromo
cd chromo

# Pull the manylinux image
docker pull quay.io/pypa/manylinux2014_x86_64

# Start container with chromo mounted
docker run --rm -d -it --name chromo -v "$(pwd)":/app quay.io/pypa/manylinux2014_x86_64

# Inside Docker
docker exec -it chromo bash
cd /app
python3.11 -m venv venv && source venv/bin/activate
pip install --prefer-binary -v -e .
```

For Apple Silicon, use the `manylinux2014_aarch64` image instead.
