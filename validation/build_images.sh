#!/usr/bin/env bash
# Build the G4beamline image and pull the BDSIM image used by the validation
# harness.  Run once before validation/compare.py.
#
# Notes
# -----
# * G4beamline 3.06 is compiled from source on top of a prebuilt Geant4 10.5
#   base image (wmoore28/geant4:10.5.0_with_source).  The source tarball is
#   fetched from the JeffersonLab/docker-g4beamline repository.
# * The base image is CentOS 7 (EOL), so the Dockerfile redirects yum to
#   vault.centos.org.  In a sandboxed network where only an HTTPS proxy is
#   reachable, build with --network=host and copy the proxy CA bundle in (see
#   the README and the PROXY build-arg in the Dockerfile).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
G4BL_DIR="$HERE/g4bl-docker"
G4BL_IMAGE="${G4BL_IMAGE:-g4beamline:3.06}"
BDSIM_IMAGE="${BDSIM_IMAGE:-bdsim/ubuntu24-g4.11.3-bdsim:develop}"
TARBALL_URL="https://raw.githubusercontent.com/JeffersonLab/docker-g4beamline/master/tarballs/G4beamline-3.06-source.tgz"

mkdir -p "$G4BL_DIR/tarballs"
if [ ! -f "$G4BL_DIR/tarballs/G4beamline-3.06-source.tgz" ]; then
    echo ">> Downloading G4beamline 3.06 source tarball ..."
    curl -fSL -o "$G4BL_DIR/tarballs/G4beamline-3.06-source.tgz" "$TARBALL_URL"
fi

# Provide the agent-proxy CA bundle if present (needed for yum->vault via proxy).
if [ -f /root/.ccr/ca-bundle.crt ]; then
    cp /root/.ccr/ca-bundle.crt "$G4BL_DIR/ca-bundle.crt"
else
    # Empty placeholder so the COPY in the Dockerfile succeeds off-sandbox.
    : > "$G4BL_DIR/ca-bundle.crt"
fi

BUILD_ARGS=()
NET_ARG=()
if [ -n "${HTTPS_PROXY:-}" ]; then
    BUILD_ARGS+=(--build-arg "PROXY=${HTTPS_PROXY}")
    NET_ARG+=(--network=host)
fi

echo ">> Building $G4BL_IMAGE ..."
docker build "${NET_ARG[@]}" "${BUILD_ARGS[@]}" -t "$G4BL_IMAGE" "$G4BL_DIR"

echo ">> Pulling $BDSIM_IMAGE ..."
docker pull "$BDSIM_IMAGE"

echo ">> Done.  Images:"
docker images | grep -E "g4beamline|bdsim" || true
