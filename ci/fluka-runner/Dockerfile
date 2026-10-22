# Self-hosted GitHub Actions runner for chromo's FLUKA tests.
#
# FLUKA is not part of the image (no redistribution rights). It is
# bind-mounted read-only at /opt/fluka (= $FLUPRO) when the container runs.
# The toolchain is the one the prebuilt FLUKA 2025.1 libraries
# (gfortran 10.3 / glibc 2.32) are known to link with: EL9 (glibc 2.34) and
# gcc-toolset-14 (gfortran 14.2).
FROM almalinux:9.8

ARG RUNNER_VERSION=2.338.0
ARG RUNNER_SHA256=af4b794c1bc41d73d40535e3fe092a39f9679cd8d965954c2aca25a05ca41d32
# FLUKA ships data files with mode 600: use the uid that owns the mounted
# FLUKA copy, or make the copy group-readable (the entrypoint adds the group
# that owns $FLUPRO)
ARG RUNNER_UID=1001

RUN dnf -y install \
        gcc-toolset-14-gcc gcc-toolset-14-gcc-c++ gcc-toolset-14-gcc-gfortran \
        git make tar gzip xz jq util-linux findutils which procps-ng \
    && dnf clean all

RUN useradd --create-home --uid "${RUNNER_UID}" runner \
    && mkdir -p /home/runner/actions-runner /home/runner/.cache
WORKDIR /home/runner/actions-runner
RUN curl -fsSL -o runner.tgz \
        "https://github.com/actions/runner/releases/download/v${RUNNER_VERSION}/actions-runner-linux-x64-${RUNNER_VERSION}.tar.gz" \
    && echo "${RUNNER_SHA256}  runner.tgz" | sha256sum -c - \
    && tar xzf runner.tgz && rm runner.tgz \
    && ./bin/installdependencies.sh && dnf clean all \
    && chown -R runner:runner /home/runner

COPY entrypoint.sh job-started.sh /usr/local/bin/

ENV PATH=/opt/rh/gcc-toolset-14/root/usr/bin:$PATH \
    LANG=C.UTF-8 \
    FLUPRO=/opt/fluka \
    ACTIONS_RUNNER_HOOK_JOB_STARTED=/usr/local/bin/job-started.sh

# starts as root only to read the GitHub token; the runner runs as `runner`
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
