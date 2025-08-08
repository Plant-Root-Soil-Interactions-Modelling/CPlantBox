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
    cmake . && make -j$(nproc)
    # Create ABI-suffixed alias for the built extension if needed
    python3 - <<PY
import glob, os, shutil
try:
    from importlib import machinery as m
    suffixes = getattr(m, "EXTENSION_SUFFIXES", [".so"])
except Exception:
    suffixes = [".so"]
expected = "plantbox" + (suffixes[0] if suffixes else ".so")
files = sorted(glob.glob("plantbox*.so"))
print("Discovered:", files)
if files and not os.path.exists(expected):
    src = files[0]
    try:
        os.link(src, expected)
        print("Hard-linked", src, "->", expected)
    except Exception:
        shutil.copy2(src, expected)
        print("Copied", src, "->", expected)
print("Final:", expected, os.path.exists(expected))
PY
    # Curated headless smoke tests (known to be stable in container)
    cd test
    OMP_NUM_THREADS=1 python3 -m pytest -q \
      test_organ.py \
      test_organparameter.py \
      test_rootparameter.py \
      test_leafparameter.py \
      test_leaf.py \
      test_hydraulic_parameters.py
    cd ..
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


