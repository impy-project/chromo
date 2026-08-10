#!/usr/bin/env bash
#
# Pre-flight check for the macOS wheel toolchain, for use before cibuildwheel.
#
# Asserts that gfortran compiles and links, then derives
# MACOSX_DEPLOYMENT_TARGET from the gfortran runtime dylibs and exports it via
# $GITHUB_ENV.  Homebrew builds those dylibs for the runner's own OS release, so
# a target hardcoded below theirs makes delocate-wheel reject the finished wheel.
#
# Optional inputs:
#   FC                        Fortran compiler to probe (default: gfortran).
#   EXPECTED_MACOS_TARGET     Minimum macOS the released wheels should support.
#                             Warns if the toolchain forces something newer.

set -euo pipefail

if [ "$(uname -s)" != "Darwin" ]; then
  echo "Not macOS ($(uname -s)) -- nothing to check."
  exit 0
fi

FC=${FC:-gfortran}
EXPECTED_MACOS_TARGET=${EXPECTED_MACOS_TARGET:-}

fail() {
  echo "::error title=macOS wheel pre-flight::$*"
  exit 1
}

# --- 1. is the Fortran compiler actually usable? ------------------------------

fc_path=$(command -v "$FC" || true)
if [ -z "$fc_path" ]; then
  fail "Fortran compiler '$FC' is not on PATH. Check the setup-fortran step."
fi
if [ ! -x "$fc_path" ]; then
  # A dangling symlink resolves as a path but is not executable.
  fail "'$FC' resolves to $fc_path, which is not executable (dangling symlink?). \
Homebrew prefix here is $(brew --prefix 2>/dev/null || echo unknown)."
fi
echo "$FC -> $fc_path"
"$FC" --version | head -1

# meson compiles and links, so probe that, not just the binary's existence.
probe_dir=$(mktemp -d)
trap 'rm -rf "$probe_dir"' EXIT
cat >"$probe_dir/probe.f90" <<'EOF'
subroutine probe(x)
  implicit none
  double precision, intent(inout) :: x
  x = sqrt(x) + 1.0d0
end subroutine probe
EOF
"$FC" -shared -fPIC -o "$probe_dir/libprobe.dylib" "$probe_dir/probe.f90" \
  || fail "'$FC' cannot compile and link a trivial Fortran shared library."
echo "Fortran compile + link OK"

# --- 2. what deployment target do the bundled runtime dylibs require? ---------

# vtool reports the LC_BUILD_VERSION / LC_VERSION_MIN_MACOSX of a Mach-O file.
minos_of() {
  local lib=$1
  vtool -show-build "$lib" 2>/dev/null |
    awk '$1 == "minos" || $1 == "version" { print $2; exit }'
}

# Everything delocate copies into <pkg>/.dylibs comes from the gcc keg.
runtime_libs=(libgfortran.dylib libquadmath.dylib libstdc++.dylib libgcc_s.1.1.dylib)

required=""
for name in "${runtime_libs[@]}"; do
  lib=$("$FC" --print-file-name="$name")
  # --print-file-name echoes the input unchanged when it finds nothing.
  [ "$lib" != "$name" ] && [ -e "$lib" ] || continue
  minos=$(minos_of "$lib")
  [ -n "$minos" ] || continue
  echo "$name ($lib) has a minimum target of $minos"
  if [ -z "$required" ] ||
     [ "$(printf '%s\n%s\n' "$required" "$minos" | sort -V | tail -1)" = "$minos" ]; then
    required=$minos
  fi
done

if [ -z "$required" ]; then
  required="$(sw_vers -productVersion | cut -d. -f1).0"
  echo "::warning title=macOS wheel pre-flight::Could not read a minimum target \
from the gfortran runtime dylibs; falling back to the runner OS ($required)."
fi

echo "Deployment target required by the toolchain: $required"

# delocate would reject the wheel over this an hour from now.
if [ -n "${MACOSX_DEPLOYMENT_TARGET:-}" ] &&
   [ "$MACOSX_DEPLOYMENT_TARGET" != "$required" ] &&
   [ "$(printf '%s\n%s\n' "$MACOSX_DEPLOYMENT_TARGET" "$required" | sort -V | tail -1)" = "$required" ]; then
  fail "MACOSX_DEPLOYMENT_TARGET=$MACOSX_DEPLOYMENT_TARGET is older than the \
$required required by the gfortran runtime dylibs on this runner; \
delocate-wheel would reject the wheel."
fi

if [ -n "$EXPECTED_MACOS_TARGET" ] && [ "$EXPECTED_MACOS_TARGET" != "$required" ]; then
  echo "::warning title=macOS support floor moved::Wheels built here will \
require macOS $required, not the expected $EXPECTED_MACOS_TARGET. The runner \
image probably moved to a newer macOS; pin the runner label or update \
EXPECTED_MACOS_TARGET deliberately."
fi

if [ -n "${GITHUB_ENV:-}" ]; then
  echo "MACOSX_DEPLOYMENT_TARGET=$required" >>"$GITHUB_ENV"
fi
