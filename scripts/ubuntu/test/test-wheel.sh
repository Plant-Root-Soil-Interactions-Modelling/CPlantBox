#!/usr/bin/env bash
set -euo pipefail

# Test an existing Ubuntu-built wheel by installing it into a venv inside a
# dedicated Ubuntu test image, then running the golden smoke.
# Usage: scripts/ubuntu/test/test-wheel.sh [linux/amd64|linux/arm64] [WHEEL_PATH(optional)]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$REPO_ROOT"

source "$REPO_ROOT/scripts/common/env.sh"
source "$REPO_ROOT/scripts/common/logging.sh"
source "$REPO_ROOT/scripts/common/docker.sh"

PLATFORM="${1:-linux/amd64}"
ARCH_SHORT="${PLATFORM#linux/}"
IMAGE_NAME="cplantbox-ubuntu-wheel-test-${ARCH_SHORT}"
if [[ "$ARCH_SHORT" == "arm64" ]]; then
  DOCKERFILE_PATH="scripts/ubuntu/test/Dockerfile_arm64"
else
  DOCKERFILE_PATH="scripts/ubuntu/test/Dockerfile_amd_64"
fi

# Find wheel if not provided
WHEEL_ARG="${2:-}"
if [[ -z "${WHEEL_ARG}" ]]; then
  WHEEL_DIR="wheelhouse/linux/ubuntu/${ARCH_SHORT}"
  if [[ ! -d "$WHEEL_DIR" ]]; then
    die "No wheelhouse dir found at $WHEEL_DIR. Run build first."
  fi
  WHEEL_ARG=$(ls -1t "$WHEEL_DIR"/*.whl | head -n1)
fi

if [[ ! -f "$WHEEL_ARG" ]]; then
  die "Wheel not found: $WHEEL_ARG"
fi

log_info "[ubuntu-test] Ensuring buildx builder..."
ensure_buildx_builder cplantbox-builder

log_info "[ubuntu-test] Building image ${IMAGE_NAME} from ${DOCKERFILE_PATH} for ${PLATFORM}..."
docker_build_image "$IMAGE_NAME" "$DOCKERFILE_PATH" "$PLATFORM" .

# Compute repo-relative wheel path so it exists inside the container mount
WHEEL_ABS=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$WHEEL_ARG")
REL_WHL="${WHEEL_ABS#${REPO_ROOT}/}"

log_info "[ubuntu-test] Running golden smoke in container venv..."
docker_run_project "$IMAGE_NAME" "$PLATFORM" /src "bash scripts/ubuntu/test/run_golden_in_container.sh '$REL_WHL'"

log_info "[ubuntu-test] Completed."


