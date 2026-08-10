#!/usr/bin/env bash
#
# Fast pre-flight check for the macOS wheel toolchain.
#
# Two failure modes have taken the Build workflow down on main without any PR
# ever seeing them.  Both are cheap to detect up front, but expensive to hit
# inside cibuildwheel, where they only surface after the Fortran models have
# been compiled:
#
#   1. gfortran is not usable.  setup-fortran <= v1.8 symlinked
#      /opt/homebrew/bin/gfortran-<v> into /usr/local/bin unconditionally, so
#      on the Intel runners (Homebrew prefix /usr/local) the symlink dangles.
#      meson then fails with "Unknown compiler(s): [['gfortran']]".
#
#   2. MACOSX_DEPLOYMENT_TARGET is older than the gfortran runtime dylibs that
#      delocate bundles into the wheel.  Homebrew builds those for the runner's
#      own OS release, so any hardcoded target breaks the moment the image moves
#      (macos-latest -> macOS 26 gives dylibs with a minimum target of 16.0).
#      delocate-wheel rejects the wheel with
#      "Library dependencies do not satisfy target MacOS version".
#
# On success the deployment target derived from the toolchain itself is written
# to $GITHUB_ENV, so the wheel job has a single source of truth instead of a
# hardcoded number that has to be remembered on every image bump.
#
# Optional inputs:
#   FC                        Fortran compiler to probe (default: gfortran).
#   EXPECTED_MACOS_TARGET     Deployment target the released wheels are meant to
#                             support.  If the toolchain forces something newer,
#                             emit a warning so a silently raised support floor
#                             does not go unnoticed.

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
Homebrew prefix here is $(brew --prefix 2>/dev/null || echo unknown); \
setup-fortran must derive the prefix rather than hardcode /opt/homebrew."
fi
echo "$FC -> $fc_path"
"$FC" --version | head -1

# Compiling and linking is what meson actually does, so prove it works rather
# than trusting that the binary exists.
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

# A target that was already pinned lower than the runtime dylibs is exactly the
# state that makes delocate fail, so say so here instead of an hour from now.
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
