#!/usr/bin/env bash
set -euo pipefail

# Audit the built Linux wheel for external shared library dependencies.
# This uses the Ubuntu test container; for full manylinux compliance, we may
# later switch to a manylinux docker image and use `auditwheel repair`.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

IMAGE_NAME="cplantbox-ubuntu-test-env"
DOCKERFILE_PATH="docker/Dockerfile.ubuntu-test-env"
PLATFORM="linux/amd64"

if ! docker buildx inspect cplantbox-builder >/dev/null 2>&1; then
  docker buildx create --use --name cplantbox-builder >/dev/null
fi

echo "[audit] Ensuring Ubuntu build image..."
docker buildx build --platform "$PLATFORM" -f "$DOCKERFILE_PATH" -t "$IMAGE_NAME" --load .

echo "[audit] Building wheel (fresh) inside container..."
docker run --rm -t \
  --platform "$PLATFORM" \
  -v "$(pwd)":/src -w /src \
  "$IMAGE_NAME" bash -lc '
    set -euo pipefail
    python3 -m pip install -U pip build auditwheel
    rm -rf dist/ _skbuild/
    python3 -m build --wheel
    echo "Built wheel:" && ls -lh dist/*.whl
    echo
    echo "[audit] auditwheel show:" && echo
    auditwheel show dist/*.whl || true
  '

echo "[audit] Done."


