#!/usr/bin/env bash
set -euo pipefail

# Build manylinux2014 wheels inside the official pypa container and auditwheel-repair them.
# Starts with cp310 to validate the flow; can extend to 3.9/3.11/3.12.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

PLATFORM="linux/amd64"
IMAGE="quay.io/pypa/manylinux2014_x86_64"
# Note: We will autodetect available Python versions inside the container.

echo "[manylinux] Pulling $IMAGE..."
docker pull "$IMAGE" >/dev/null

echo "[manylinux] Building wheels (autodetecting /opt/python cp39–cp312)"
docker run --rm -t \
  --platform "$PLATFORM" \
  -v "$(pwd)":/project -w /project \
  "$IMAGE" bash -lc '
    set -euo pipefail
    mkdir -p wheelhouse
    for PYBIN in /opt/python/cp3{9,10,11,12}-*/bin; do
      [ -d "$PYBIN" ] || continue
      PY=$(basename $(dirname "$PYBIN"))
      cd /project
      "$PYBIN/python" -m pip install -U pip build scikit-build-core cmake ninja auditwheel
      rm -rf dist/ _skbuild/
      export CMAKE_ARGS="-DPython_EXECUTABLE=$PYBIN/python -DPython_ROOT_DIR=${PYBIN%/bin}"
      "$PYBIN/python" -m build --wheel
      auditwheel show dist/*.whl || true
      auditwheel repair -w wheelhouse dist/*.whl || true
      echo "[manylinux] Testing repaired wheel for $PY..."
      "$PYBIN/python" -m pip install -U virtualenv >/dev/null
      "$PYBIN/virtualenv" -p "$PYBIN/python" /tmp/venv_test
      /tmp/venv_test/bin/python -m pip install -U pip >/dev/null
      WHL=$(ls -1 wheelhouse/*-${PY}-*.whl | head -n1)
      /tmp/venv_test/bin/python -m pip install "$WHL" >/dev/null
      cd /tmp && /tmp/venv_test/bin/python - <<PY
import os
import plantbox as pb
rs = pb.Plant()
rootsys = os.path.join(pb.data_path(), "structural", "rootsystem", "Anagallis_femina_Leitner_2010.xml")
rs.readParameters(rootsys, verbose=False)
rs.initialize(False)
rs.simulate(1.0, False)
ext_path = pb._ext.__file__ if hasattr(pb, "_ext") else getattr(pb, "__file__", "?")
try:
    ext_size = os.path.getsize(pb._ext.__file__)
except Exception:
    ext_size = -1
print("Extension:", ext_path)
print("Extension size(bytes):", ext_size)
print("manylinux wheel import/simulate OK")
PY
      cd /project
    done
     echo "Wheelhouse contents (on /project):" && ls -lh /project/wheelhouse/ || true
  '

echo "[manylinux] Done. Artifacts in wheelhouse/"


