#!/usr/bin/env bash
set -euo pipefail

# Simple helper to build the Ubuntu (amd64) testing image with buildx and run tests
# against the current working tree (mounted into /src).

IMAGE_NAME="cplantbox-ubuntu-test-env"
DOCKERFILE_PATH="docker/Dockerfile.ubuntu-test-env"
PLATFORM="linux/amd64"

# Ensure buildx is set up
if ! docker buildx inspect cplantbox-builder >/dev/null 2>&1; then
  docker buildx create --use --name cplantbox-builder >/dev/null
fi

echo "Building $IMAGE_NAME from $DOCKERFILE_PATH for $PLATFORM..."
docker buildx build --platform "$PLATFORM" -f "$DOCKERFILE_PATH" -t "$IMAGE_NAME" --load .

echo "Running tests inside container..."
docker run --rm -t \
  --platform "$PLATFORM" \
  -v "$(pwd)":/src -w /src \
  "$IMAGE_NAME" bash -lc '
    set -euo pipefail
    export PYTHONPATH="/src:/src/src:${PYTHONPATH:-}"
    git submodule update --init --recursive

    echo "Attempting wheel build and wheel-based smoke..."
    if python3 -m pip -q install build && rm -rf dist/ _skbuild/ && SETUPTOOLS_SCM_PRETEND_VERSION="$(git describe --tags --abbrev=0 2>/dev/null || echo 2.1.0)" python3 -m build --wheel >/dev/null 2>&1; then
      echo "Wheel built; testing installed package..."
      python3 -m venv /tmp/venv && . /tmp/venv/bin/activate
      python3 -m pip -q install --upgrade pip
      WHEEL_PATH=$(ls -1t dist/*.whl | head -n1)
      python3 -m pip -q install "$WHEEL_PATH"
      python3 - <<PY
import plantbox as pb
rs = pb.Plant()
rs.readParameters("modelparameter/structural/rootsystem/Anagallis_femina_Leitner_2010.xml", verbose=False)
rs.initialize(False)
rs.simulate(5.0, False)
print("Wheel import/simulate OK")
PY
      deactivate
    else
      echo "Wheel build failed; falling back to source-tree smoke..."
      rm -f CMakeCache.txt && rm -rf CMakeFiles && find . -maxdepth 1 -name "*.so" -delete
      cmake . && make -j$(nproc)
      cd test
      OMP_NUM_THREADS=1 python3 -m pytest -q \
        test_organ.py \
        test_organparameter.py \
        test_rootparameter.py \
        test_leafparameter.py \
        test_leaf.py \
        test_hydraulic_parameters.py
      cd ..
    fi

    # Headless smoke scenario: build a small plant and simulate briefly (no VTK)
    python3 - <<PY
import plantbox as pb
rs = pb.Plant()
rs.readParameters("modelparameter/structural/rootsystem/Anagallis_femina_Leitner_2010.xml", verbose=False)
rs.initialize(False)
rs.simulate(5.0, False)
rs.write("results/smoke_example.vtp")
print("Smoke simulation OK; wrote results/smoke_example.vtp")
PY
  '

echo "All tests completed."


