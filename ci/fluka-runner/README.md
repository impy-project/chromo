# FLUKA CI runner

A self-hosted GitHub Actions runner that builds chromo against a local FLUKA
installation and runs the FLUKA tests (`.github/workflows/fluka.yml`). FLUKA
has no redistribution rights, so it is not in the image or on GitHub: it is
mounted read-only from the runner host, and nothing built from it leaves the
runner (the workflow uploads no artifacts and uses no caches).

## How it works

- Each job runs in a fresh container (`docker compose run --rm`, restarted
  by systemd). `entrypoint.sh` registers an ephemeral runner with the label
  `fluka`, runs one job, and deregisters.
- The toolchain matches the one the prebuilt FLUKA 2025.1 libraries
  (gfortran 10.3, glibc 2.32) are known to link with: AlmaLinux 9.8
  (glibc 2.34) and gcc-toolset-14 (gfortran 14.2).
- The entrypoint starts as root only to read the GitHub token, which jobs
  cannot read, and runs the runner and the jobs as the unprivileged user
  `runner`.
- `job-started.sh` runs before every job and refuses anything except
  `fluka.yml` of `impy-project/chromo` from a branch, started by `push` or
  `workflow_dispatch`. The repository is public: a fork pull request can add
  its own workflow with `runs-on: fluka`, and on GitHub Free there are no
  runner groups to restrict which workflows may use a runner, so the runner
  enforces it. Pull request branches are tested with
  `gh workflow run fluka.yml --ref <branch>`, which needs write access.

## Setup on the runner host

1. Docker with the compose plugin, and a FLUKA 2025.1 installation
   (`scripts/install_fluka.sh`; the patched `libflukahp.a` for specific
   isotopes and free-neutron targets).
2. FLUKA ships data files with mode 600. Either build the image with
   `RUNNER_UID` set to the owner of the FLUKA copy, or make the copy readable
   by its group (`chmod -R g+rX`); the entrypoint adds the group that owns
   `$FLUPRO` to `runner`. The entrypoint stops with a message if any file is
   unreadable, since meson would otherwise skip the FLUKA build.
3. A token that can create runner registration tokens for the repository:
   a fine-grained personal access token for `impy-project/chromo` with
   **Administration: read and write**. Store it as `gh_token` next to
   `compose.yaml`, mode 600 (it is git-ignored).
4. Build and start:

   ```bash
   sudo mkdir -p /srv/chromo-fluka-runner && sudo cp -r ci/fluka-runner/. /srv/chromo-fluka-runner/
   cd /srv/chromo-fluka-runner
   RUNNER_UID=$(stat -c %u /srv/fluka/2025.1) docker compose build
   # edit WorkingDirectory and FLUKA_DIR in the unit, then
   sudo cp fluka-runner.service /etc/systemd/system/
   sudo systemctl enable --now fluka-runner
   ```

5. In the repository settings, keep "Require approval for all external
   contributors" for fork pull request workflows.

## Maintenance

- The runner runs with `--disableupdate`. GitHub stops sending jobs to
  runners that are more than 30 days behind the latest release, so bump
  `RUNNER_VERSION` and `RUNNER_SHA256` in the `Dockerfile` (the checksum is
  in the release notes of `actions/runner`) and rebuild.
- `test_generators` has no committed FLUKA references; on a fresh checkout
  those cases write new references and report xfail.
