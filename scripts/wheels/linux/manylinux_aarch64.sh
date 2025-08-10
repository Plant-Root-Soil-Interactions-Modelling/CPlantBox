#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$REPO_ROOT"

source "$REPO_ROOT/scripts/common/env.sh"
source "$REPO_ROOT/scripts/common/logging.sh"

IMAGE="quay.io/pypa/manylinux2014_aarch64"
PLATFORM="linux/arm64"
OUT_DIR="wheelhouse/linux/manylinux2014_aarch64"
mkdir -p "$OUT_DIR"

log_info "[manylinux aarch64] Pulling $IMAGE..."
docker pull "$IMAGE" >/dev/null

log_info "[manylinux aarch64] Building + repairing wheels (cp39-cp312)"
docker run --rm -t \
  --platform "$PLATFORM" \
  -v "$(pwd)":/project -w /project \
  "$IMAGE" bash -lc $'\
    set -euo pipefail\n\
    mkdir -p /project/'"$OUT_DIR"$'\n\
    for PYBIN in /opt/python/cp3{9,10,11,12}-*/bin; do\n\
      [ -d "$PYBIN" ] || continue\n\
      PY=$(basename $(dirname "$PYBIN"))\n\
      echo "[manylinux aarch64] Building for $PY"\n\
      cd /project\n\
      "$PYBIN/python" -m pip install -U pip build scikit-build-core cmake ninja auditwheel >/dev/null\n\
      rm -rf dist/ _skbuild/\n\
      export CMAKE_ARGS="-DPython_EXECUTABLE=$PYBIN/python -DPython_ROOT_DIR=${PYBIN%/bin}"\n\
      "$PYBIN/python" -m build --wheel\n\
      auditwheel show dist/*.whl || true\n\
      auditwheel repair -w /project/'"$OUT_DIR"$' dist/*.whl || true\n\
      echo "[manylinux aarch64] Smoke test for $PY"\n\
      "$PYBIN/python" -m pip install -U virtualenv >/dev/null\n\
      "$PYBIN/virtualenv" -p "$PYBIN/python" /tmp/venv_test\n\
      /tmp/venv_test/bin/python -m pip install -U pip >/dev/null\n\
      WHL=$(ls -1 /project/'"$OUT_DIR"$'/*-${PY}-*.whl | head -n1)\n\
      /tmp/venv_test/bin/python -m pip install "$WHL" >/dev/null\n\
      cd /tmp && /tmp/venv_test/bin/python - <<PY\nimport os, plantbox as pb\nprint("Imported plantbox", getattr(pb, "__version__", "?"))\nrootsys = os.path.join(pb.data_path(), "structural", "rootsystem", "Anagallis_femina_Leitner_2010.xml")\nprint("Param exists:", os.path.exists(rootsys))\nrs = pb.Plant()\nrs.readParameters(rootsys, verbose=False)\nrs.initialize(False)\nrs.simulate(0.5, False)\nprint("manylinux aarch64 wheel import/simulate OK")\nPY\n\
      rm -rf /tmp/venv_test\n\
      cd /project\n\
    done\n\
    echo "[manylinux aarch64] Artifacts:" && ls -lh /project/'"$OUT_DIR"$' || true\n'

log_info "[manylinux aarch64] Done. Artifacts in $OUT_DIR"


