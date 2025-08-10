#!/usr/bin/env bash
set -euo pipefail

# Build wheel inside the Ubuntu amd64 container and smoke-test it from a fresh venv

IMAGE_NAME="cplantbox-ubuntu-test-env"
DOCKERFILE_PATH="docker/Dockerfile.ubuntu-test-env"
PLATFORM="linux/amd64"

if ! docker buildx inspect cplantbox-builder >/dev/null 2>&1; then
  docker buildx create --use --name cplantbox-builder >/dev/null
fi

echo "Building $IMAGE_NAME from $DOCKERFILE_PATH for $PLATFORM..."
docker buildx build --platform "$PLATFORM" -f "$DOCKERFILE_PATH" -t "$IMAGE_NAME" --load .

echo "Building wheel and testing install inside container..."
docker run --rm -t \
  --platform "$PLATFORM" \
  -v "$(pwd)":/src -w /src \
  "$IMAGE_NAME" bash -lc '
    set -euo pipefail
    python3 -m pip install -U pip build
    # Determine version from git tag if available; fallback to 2.1.0
    VERSION="$(git describe --tags --abbrev=0 2>/dev/null || echo 2.1.0)"
    export SETUPTOOLS_SCM_PRETEND_VERSION="$VERSION"
    rm -rf dist/ _skbuild/
    python3 -m build --wheel
    echo "Built wheel:" && ls -lh dist/*.whl
    python3 -m venv /tmp/venv
    /tmp/venv/bin/python -m pip install -U pip
    WHEEL_PATH=$(ls -1t dist/*.whl | head -n1)
    /tmp/venv/bin/python -m pip install "$WHEEL_PATH"
    cd /tmp
    /tmp/venv/bin/python - <<PY
import os
import plantbox as pb
rs = pb.Plant()
rootsys = os.path.join(pb.data_path(), "structural", "rootsystem", "Anagallis_femina_Leitner_2010.xml")
rs.readParameters(rootsys, verbose=False)
rs.initialize(False)
rs.simulate(2.0, False)
ext_path = pb._ext.__file__ if hasattr(pb, "_ext") else getattr(pb, "__file__", "?")
try:
    ext_size = os.path.getsize(pb._ext.__file__)
except Exception:
    ext_size = -1
print("Extension:", ext_path)
print("Extension size(bytes):", ext_size)
print("Wheel import+simulate OK with packaged data")
PY
  '

echo "Wheel install smoke OK."


