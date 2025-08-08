#!/usr/bin/env bash
set -euo pipefail

# Build manylinux2014 wheels inside the official pypa container and auditwheel-repair them.
# Starts with cp310 to validate the flow; can extend to 3.9/3.11/3.12.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

PLATFORM="linux/amd64"
IMAGE="quay.io/pypa/manylinux2014_x86_64"
PY_TAGS=(cp310-cp310)

echo "[manylinux] Pulling $IMAGE..."
docker pull "$IMAGE" >/dev/null

echo "[manylinux] Building wheels for: ${PY_TAGS[*]}"
docker run --rm -t \
  --platform "$PLATFORM" \
  -v "$(pwd)":/project -w /project \
  "$IMAGE" bash -lc '
    set -euo pipefail
    mkdir -p wheelhouse
    for PY in ${PY_TAGS[@]:-cp310-cp310}; do
      PYBIN="/opt/python/${PY}/bin"
      if [ ! -x "$PYBIN/python" ]; then
        echo "Skipping $PY (not present)" >&2
        continue
      fi
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
      WHL=$(ls -1 wheelhouse/*.whl | head -n1)
      /tmp/venv_test/bin/python -m pip install "$WHL" >/dev/null
      /tmp/venv_test/bin/python - <<PY
import os
import plantbox as pb
rs = pb.Plant()
rootsys = os.path.join(pb.data_path(), "structural", "rootsystem", "Anagallis_femina_Leitner_2010.xml")
rs.readParameters(rootsys, verbose=False)
rs.initialize(False)
rs.simulate(1.0, False)
print("manylinux wheel import/simulate OK")
PY
    done
    echo "Wheelhouse contents:" && ls -lh wheelhouse/
  '

echo "[manylinux] Done. Artifacts in wheelhouse/"


