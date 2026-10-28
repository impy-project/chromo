#!/usr/bin/env bash
# ACTIONS_RUNNER_HOOK_JOB_STARTED: runs before any step of a job; a non-zero
# exit fails the job.
#
# The repository is public, and a pull request from a fork can add its own
# workflow with `runs-on: fluka`, which would expose the mounted FLUKA
# installation. Only accept this repository's fluka.yml from a branch, started
# by a push or by a user with write access (workflow_dispatch). Unset
# variables fail the hook, so the check fails closed.
set -euo pipefail

refuse() {
    echo "FLUKA runner: refusing job ($1): repository=${GITHUB_REPOSITORY:-?}" \
        "event=${GITHUB_EVENT_NAME:-?} workflow_ref=${GITHUB_WORKFLOW_REF:-?}" >&2
    exit 1
}

[[ "$GITHUB_REPOSITORY" == "$GH_REPOSITORY" ]] || refuse "repository"
case "$GITHUB_EVENT_NAME" in
    push | workflow_dispatch) ;;
    *) refuse "event" ;;
esac
[[ "$GITHUB_WORKFLOW_REF" == "$GH_REPOSITORY/.github/workflows/fluka.yml@refs/heads/"* ]] \
    || refuse "workflow"

echo "FLUKA runner: accepted $GITHUB_WORKFLOW_REF ($GITHUB_EVENT_NAME by ${GITHUB_ACTOR:-?})"
