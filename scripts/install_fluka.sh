#!/usr/bin/env bash
# Install FLUKA 2025.1 into $FLUPRO (default: $HOME/devel/FLUKA).
# Idempotent: if $FLUPRO/libflukahp.a exists, only verifies the install.

set -euo pipefail

FLUPRO="${FLUPRO:-$HOME/devel/FLUKA}"
ARCHIVE_DIR="${FLUKA_ARCHIVE_DIR:-$HOME/devel/FLUKA-dev}"
CODE_TAR="$ARCHIVE_DIR/fluka2025.1-linux-gfor64bit-10.3-glibc2.32-AA.tar.gz"
DATA_TAR="$ARCHIVE_DIR/fluka2025.1-data.tar.gz"

REQUIRED=(
  "libflukahp.a"
  "libdpmmvax.a"
  "librqmdmvax.a"
  "latestRQMD/librqmd.a"
  "interface/dpmvers"
)

echo "FLUPRO=$FLUPRO"
mkdir -p "$FLUPRO"

if [[ -f "$FLUPRO/libflukahp.a" ]]; then
  echo "libflukahp.a present, skipping untar/build; verifying only."
else
  for tar in "$CODE_TAR" "$DATA_TAR"; do
    [[ -f "$tar" ]] || { echo "ERROR: archive not found: $tar" >&2; exit 1; }
  done
  echo "Extracting archives into $FLUPRO..."
  ( cd "$FLUPRO" && tar xzf "$CODE_TAR" && tar xzf "$DATA_TAR" )
  echo "Running FLUKA top-level make..."
  ( cd "$FLUPRO" && make )
fi

for f in "${REQUIRED[@]}"; do
  if [[ ! -f "$FLUPRO/$f" ]]; then
    echo "ERROR: missing expected file: $FLUPRO/$f" >&2
    exit 1
  fi
done
ls "$FLUPRO"/interface/libdpmjet*.a >/dev/null 2>&1 || {
  echo "ERROR: no libdpmjet*.a under $FLUPRO/interface/" >&2; exit 1;
}

echo ""
echo "FLUKA install OK at $FLUPRO"
echo "Add to your shell rc:"
echo "  export FLUPRO=$FLUPRO"
