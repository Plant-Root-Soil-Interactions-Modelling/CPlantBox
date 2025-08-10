#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

source "$REPO_ROOT/scripts/common/env.sh"
source "$REPO_ROOT/scripts/common/logging.sh"
source "$REPO_ROOT/scripts/common/docker.sh"

IMAGE_NAME="cplantbox-ubuntu-test-env"
DOCKERFILE_PATH="docker/Dockerfile.ubuntu-test-env"
PLATFORM="linux/amd64"

log_info "Ensuring buildx builder..."
ensure_buildx_builder cplantbox-builder

log_info "Building $IMAGE_NAME from $DOCKERFILE_PATH for $PLATFORM..."
docker_build_image "$IMAGE_NAME" "$DOCKERFILE_PATH" "$PLATFORM" .

log_info "Image ready: $IMAGE_NAME ($PLATFORM)"


