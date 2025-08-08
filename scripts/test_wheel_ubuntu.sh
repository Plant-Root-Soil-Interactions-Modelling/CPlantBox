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
    rm -rf dist/ _skbuild/
    python3 -m build --wheel
    echo "Built wheel:" && ls -lh dist/*.whl
    python3 -m venv /tmp/venv
    /tmp/venv/bin/python -m pip install -U pip
    WHEEL_PATH=$(ls -1t dist/*.whl | head -n1)
    /tmp/venv/bin/python -m pip install "$WHEEL_PATH"
    cd /tmp
    /tmp/venv/bin/python - <<PY
import plantbox as pb
rs = pb.Plant()
rs.readParameters("/src/modelparameter/structural/rootsystem/Anagallis_femina_Leitner_2010.xml", verbose=False)
rs.initialize(False)
rs.simulate(2.0, False)
print("Wheel import+simulate OK")
PY
  '

echo "Wheel install smoke OK."


