#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

source "$REPO_ROOT/scripts/common/env.sh"
source "$REPO_ROOT/scripts/common/logging.sh"
source "$REPO_ROOT/scripts/common/docker.sh"

IMAGE_NAME="cplantbox-ubuntu-test-env"
PLATFORM="linux/amd64"

log_info "Running wheel build + smoke inside $IMAGE_NAME ($PLATFORM)"

docker_run_project "$IMAGE_NAME" "$PLATFORM" /src $'\
  set -euo pipefail\n\
  export PYTHONPATH="/src:/src/src:${PYTHONPATH:-}"\n\
  git submodule update --init --recursive\n\
  echo "[ubuntu] Building wheel..."\n\
  python3 -m pip -q install -U pip build\n\
  rm -rf dist/ _skbuild/\n\
  python3 -m build --wheel\n\
  echo "[ubuntu] Wheel built:" && ls -lh dist/*.whl\n\
  echo "[ubuntu] Installing into temp venv and running minimal simulate..."\n\
  python3 -m venv /tmp/venv && . /tmp/venv/bin/activate\n\
  python3 -m pip -q install -U pip\n\
  WHEEL_PATH=$(ls -1t dist/*.whl | head -n1)\n\
  python3 -m pip -q install "$WHEEL_PATH"\n\
  python3 - <<PY\nimport os, plantbox as pb\nprint("Imported plantbox", getattr(pb, "__version__", "?"))\nrootsys = os.path.join(pb.data_path(), "structural", "rootsystem", "Anagallis_femina_Leitner_2010.xml")\nprint("Param exists:", os.path.exists(rootsys))\nrs = pb.Plant()\nrs.readParameters(rootsys, verbose=False)\nrs.initialize(False)\nrs.simulate(2.0, False)\nprint("Wheel import/simulate OK")\nPY\n\
  deactivate\n\
  echo "[ubuntu] Installing VTK and running golden image comparison..."\n\
  python3 -m venv /tmp/venv2 && . /tmp/venv2/bin/activate\n\
  python3 -m pip -q install -U pip pytest vtk\n\
  cp test/golden/linux/example_plant_headless.png test/golden/example_plant_headless.png\n\
  OMP_NUM_THREADS=1 python3 -m pytest -q test/test_golden_headless.py\n\
  rm -f test/golden/example_plant_headless.png\n\
  deactivate\n'

log_info "Wheel smoke finished."

# Copy built wheel(s) into a non-distribution artifact folder for convenience
OUT_DIR="wheelhouse/linux/ubuntu"
mkdir -p "$OUT_DIR"
if ls dist/*.whl >/dev/null 2>&1; then
  cp -v dist/*.whl "$OUT_DIR"/
  {
    echo "# Ubuntu-built (non-manylinux) wheels — for local testing only"
    date -u +"%Y-%m-%dT%H:%M:%SZ"
    python3 -V || true
    echo "Artifacts:"
    ls -1 "$OUT_DIR"/*.whl || true
  } > "$OUT_DIR/index.txt"
  log_info "Copied wheels to $OUT_DIR (non-portable)."
else
  log_warn "No dist/*.whl found to copy."
fi


