#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$REPO_ROOT"

source "$REPO_ROOT/scripts/common/env.sh"
source "$REPO_ROOT/scripts/common/logging.sh"

usage() {
  cat <<USAGE
Usage: $0 [x86_64|aarch64] [WHEEL_DIR|WHEEL_FILE ...]

Smoke-test repaired manylinux wheels inside the corresponding manylinux image.

Args:
  arch:           one of x86_64 (default) or aarch64
  WHEEL_DIR:      optional directory containing repaired wheels to test
  WHEEL_FILE...:  optional explicit wheel paths to test (overrides WHEEL_DIR)

If neither WHEEL_DIR nor WHEEL_FILE is given, defaults to:
  wheelhouse/linux/manylinux2014_<arch>/
USAGE
}

ARCH=${1:-x86_64}
shift || true

case "$ARCH" in
  x86_64) IMAGE="quay.io/pypa/manylinux2014_x86_64"; PLATFORM="linux/amd64"; SUBDIR="manylinux2014_x86_64" ;;
  aarch64) IMAGE="quay.io/pypa/manylinux2014_aarch64"; PLATFORM="linux/arm64"; SUBDIR="manylinux2014_aarch64" ;;
  -h|--help) usage; exit 0 ;;
  *) log_error "Unknown arch: $ARCH"; usage; exit 2 ;;
esac

WHEELS=()
if [[ $# -gt 0 ]]; then
  # Accept a directory or explicit wheel files
  if [[ $# == 1 && -d $1 ]]; then
    # Prefer newest first
    while IFS= read -r f; do [[ -n "$f" ]] && WHEELS+=("$f"); done < <(ls -1t "$1"/*.whl 2>/dev/null || true)
  else
    for f in "$@"; do WHEELS+=("$f"); done
  fi
else
  DEFAULT_DIR="wheelhouse/linux/${SUBDIR}"
  if [[ -d "$DEFAULT_DIR" ]]; then
    while IFS= read -r f; do [[ -n "$f" ]] && WHEELS+=("$f"); done < <(ls -1t "$DEFAULT_DIR"/*.whl 2>/dev/null || true)
  fi
fi

if [[ ${#WHEELS[@]} -eq 0 ]]; then
  log_warn "No wheels found to smoke. Build manylinux wheels first."
  exit 0
fi

# Limit to a single wheel by default (override with SMOKE_MAX)
SMOKE_MAX=${SMOKE_MAX:-1}
if (( ${#WHEELS[@]} > SMOKE_MAX )); then
  WHEELS=( "${WHEELS[@]:0:$SMOKE_MAX}" )
fi
log_info "[manylinux $ARCH] Will smoke ${#WHEELS[@]} wheel(s):"; for w in "${WHEELS[@]}"; do echo "  - $w"; done

log_info "[manylinux $ARCH] Pulling $IMAGE..."
docker pull "$IMAGE" >/dev/null

for WHL in "${WHEELS[@]}"; do
  if [[ ! -f "$WHL" ]]; then
    log_warn "Skipping missing: $WHL"
    continue
  fi
  ABS_WHL="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$WHL")"
  REL_WHL="$(python3 -c 'import os,sys; rp=os.path.relpath(sys.argv[1], sys.argv[2]); print(rp)' "$ABS_WHL" "$REPO_ROOT")"
  log_info "[manylinux $ARCH] Smoking wheel: $REL_WHL"

  docker run --rm -t \
    --platform "$PLATFORM" \
    -e SMOKE_GOLDEN=${SMOKE_GOLDEN:-1} \
    -v "$REPO_ROOT":/project -w /project \
    "$IMAGE" bash -lc $'\
      set -euo pipefail\n\
      WHEEL_PATH='"$REL_WHL"$'\n\
      if [[ ! -f "$WHEEL_PATH" ]]; then echo "Wheel not found: $WHEEL_PATH" >&2; exit 3; fi\n\
      FNAME=$(basename "$WHEEL_PATH")\n\
      # Extract python tag (e.g., cp310-cp310) using bash regex\n\
      if [[ "$FNAME" =~ -(cp[0-9]{2,3}-cp[0-9]{2,3})- ]]; then\n\
        PYTAG="${BASH_REMATCH[1]}"\n\
      else\n\
        echo "Cannot determine python tag for $FNAME" >&2; exit 4\n\
      fi\n\
      echo "[smoke] Python tag: $PYTAG"\n\
      # Find matching interpreter in /opt/python\n\
      MATCH_BIN=""\n\
      for p in /opt/python/${PYTAG}*/bin; do\n\
        if [[ -x "$p/python" ]]; then MATCH_BIN="$p/python"; break; fi\n\
      done\n\
      if [[ -z "$MATCH_BIN" ]]; then echo "No matching /opt/python for $PYTAG; skipping this wheel" >&2; exit 0; fi\n\
      echo "[smoke] Using interpreter: $MATCH_BIN"\n\
      # Create temp venv, install and run minimal simulate\n\
      "$MATCH_BIN" -m pip install -U pip >/dev/null\n\
      "$MATCH_BIN" -m pip install -U virtualenv >/dev/null\n\
      "$MATCH_BIN" -m virtualenv /tmp/venv_test >/dev/null\n\
      /tmp/venv_test/bin/python -m pip install -U pip >/dev/null\n\
      /tmp/venv_test/bin/python -m pip install "/project/$WHEEL_PATH" >/dev/null\n\
      cd /tmp && /tmp/venv_test/bin/python - <<PY\nimport os, plantbox as pb\nprint("Imported plantbox", getattr(pb, "__version__", "?"))\nrootsys = os.path.join(pb.data_path(), "structural", "rootsystem", "Anagallis_femina_Leitner_2010.xml")\nprint("Param exists:", os.path.exists(rootsys))\nrs = pb.Plant()\nrs.readParameters(rootsys, verbose=False)\nrs.initialize(False)\nrs.simulate(0.5, False)\nprint("manylinux wheel import/simulate OK")\nPY\n\
      if [[ "${SMOKE_GOLDEN:-1}" == "1" ]]; then\n\
        set +e\n\
        /tmp/venv_test/bin/python -m pip install -q vtk pillow numpy >/dev/null 2>&1 || true\n\
        echo "[smoke] (golden) Running headless golden test (single file)"\n\
        cp /project/test/golden/linux/example_plant_headless.png /project/test/golden/example_plant_headless.png\n\
        OMP_NUM_THREADS=1 /tmp/venv_test/bin/python /project/scripts/wheels/linux/golden_smoke.py\n\
        rm -f /project/test/golden/example_plant_headless.png\n\
        set -e\n\
      else\n\
        echo "[smoke] Golden test disabled for manylinux smoke (set SMOKE_GOLDEN=1 to enable)"\n\
      fi\n\
      rm -rf /tmp/venv_test\n'
done

log_info "[manylinux $ARCH] Smoke finished."


