#!/usr/bin/env bash
# Register an ephemeral runner, run one job and deregister.
#
# Runs as root to read the token file (mode 0400, not readable by the
# `runner` user that executes jobs). The token stays in a non-exported
# shell variable of this root process.
set -euo pipefail

: "${GH_REPOSITORY:?set GH_REPOSITORY, e.g. impy-project/chromo}"
token_file=${GH_TOKEN_FILE:-/run/secrets/gh_token}
api=${GITHUB_API_URL:-https://api.github.com}
name=${RUNNER_NAME:-fluka}-$(hostname)

if [[ ! -r "$FLUPRO/libflukahp.a" ]]; then
    echo "FLUKA is not mounted at FLUPRO=$FLUPRO" >&2
    exit 1
fi
gh_token=$(<"$token_file")

# FLUKA installations are often group-readable only; let `runner` in
# through the group that owns $FLUPRO instead of loosening permissions
fluka_gid=$(stat -c %g "$FLUPRO")

runner_token() {  # registration-token | remove-token
    curl -fsS -X POST \
        -H "Authorization: Bearer $gh_token" \
        -H "Accept: application/vnd.github+json" \
        "$api/repos/$GH_REPOSITORY/actions/runners/$1" | jq -er .token
}

as_runner() {
    HOME=/home/runner setpriv --reuid=runner --regid=runner \
        --groups="$fluka_gid" "$@"
}

# meson skips the FLUKA build when a library is unreadable, and FLUKA stops
# on an unreadable data file; fail here instead
unreadable=$(as_runner find "$FLUPRO" \
    \( \( -type f ! -readable \) -o \( -type d ! -executable \) \) -print -quit)
if [[ -n "$unreadable" ]]; then
    echo "FLUKA files are not readable by the runner user, e.g. $unreadable;" \
        "set RUNNER_UID to the owner of the FLUKA copy or make it group-readable" >&2
    exit 1
fi

reg_token=$(runner_token registration-token)
cd /home/runner/actions-runner
as_runner ./config.sh --unattended \
    --url "https://github.com/$GH_REPOSITORY" \
    --token "$reg_token" \
    --name "$name" --labels "${RUNNER_LABELS:-fluka}" \
    --work /home/runner/_work \
    --ephemeral --disableupdate --replace

# An ephemeral runner deregisters itself after its job; this covers a stop
# before that.
cleanup() {
    as_runner ./config.sh remove --token "$(runner_token remove-token)" \
        >/dev/null 2>&1 || true
}
trap cleanup EXIT

as_runner ./run.sh &
child=$!
trap 'kill -TERM "$child" 2>/dev/null || true' TERM INT
wait "$child" || true
wait "$child" 2>/dev/null || true
