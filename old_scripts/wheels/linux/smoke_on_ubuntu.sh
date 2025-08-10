#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$REPO_ROOT"

source "$REPO_ROOT/scripts/common/env.sh"
source "$REPO_ROOT/scripts/common/logging.sh"

IMAGE_NAME="cplantbox-ubuntu-test-env"
PLATFORM="linux/amd64"

WHEEL_INPUT=${1:-}

# Ensure ubuntu test image is available (also builds if missing)
"$REPO_ROOT/scripts/ubuntu/build_image.sh" >/dev/null

# Detect container Python cp tag (e.g., cp310)
CP_TAG=$(docker run --rm --platform "$PLATFORM" "$IMAGE_NAME" bash -lc "python3 -c 'import sys; print(f\"cp{sys.version_info.major}{sys.version_info.minor}\")'" 2>/dev/null || true)
if [[ -z "${CP_TAG:-}" ]]; then
  die "Could not determine container Python tag"
fi
log_info "[ubuntu-smoke] Container Python tag: $CP_TAG"

# Resolve a matching wheel to test
resolve_matching_wheel() {
  local spec="$1"  # file or directory or empty
  local match_glob="*-${CP_TAG}-*.whl"
  if [[ -n "$spec" && -f "$spec" ]]; then
    local base
    base=$(basename "$spec")
    if [[ "$base" == $match_glob ]]; then
      echo "$spec"; return 0
    else
      return 2  # mismatched provided wheel
    fi
  fi
  if [[ -n "$spec" && -d "$spec" ]]; then
    ls -1t "$spec"/$match_glob 2>/dev/null | head -n1 || true
    return 0
  fi
  # default location
  local defdir="wheelhouse/linux/manylinux2014_x86_64"
  ls -1t "$defdir"/$match_glob 2>/dev/null | head -n1 || true
}

HOST_WHEEL_PATH=$(resolve_matching_wheel "$WHEEL_INPUT" || true)
if [[ -z "${HOST_WHEEL_PATH:-}" ]]; then
  if [[ -n "$WHEEL_INPUT" && -f "$WHEEL_INPUT" ]]; then
    die "Provided wheel does not match container Python ($CP_TAG): $WHEEL_INPUT"
  fi
  die "No matching wheel (*-${CP_TAG}-*.whl) found. Build or select the correct wheel."
fi

log_info "[ubuntu-smoke] Using wheel: $HOST_WHEEL_PATH"

# Compute repo-relative path for container mount
REL_WHEEL=$(python3 -c 'import os,sys; print(os.path.relpath(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])))' "$HOST_WHEEL_PATH" "$REPO_ROOT")
if [[ -z "${REL_WHEEL:-}" || ! -f "$REL_WHEEL" ]]; then
  die "Resolved relative wheel path invalid: $REL_WHEEL"
fi

# Run golden on Ubuntu container using the installed wheel from venv; avoid importing source by running in /tmp
docker run --rm -t \
  --platform "$PLATFORM" \
  -e REL_WHEEL="$REL_WHEEL" \
  -e LIBGL_ALWAYS_SOFTWARE=1 \
  -v "$REPO_ROOT":/src -w /tmp \
  "$IMAGE_NAME" bash -lc $'\
  set -euo pipefail\n\
  python3 -m venv --system-site-packages /tmp/venv && source /tmp/venv/bin/activate\n\
  python -m pip -q install -U pip\n\
  # Pin a VTK known to behave better under Xvfb in headless testing\n  python -m pip -q install "vtk==9.2.6" || python -m pip -q install vtk\n\
  python -m pip -q install "/src/$REL_WHEEL"\n\
  # Prepare golden file\n  cp /src/test/golden/linux/example_plant_headless.png /src/test/golden/example_plant_headless.png\n\
  echo "[ubuntu-smoke] Running golden script against installed wheel under Xvfb"\n\
  set +e\n\
  xvfb-run -a -s "-screen 0 1280x1024x24" bash -lc "OMP_NUM_THREADS=1 python /src/scripts/wheels/linux/golden_smoke.py"\n  status=$?\n  set -e\n\
  exit $status\n'

status=$?
if [[ $status -ne 0 ]]; then
  log_error "[ubuntu-smoke] Golden test FAILED"
  exit $status
fi

log_info "[ubuntu-smoke] Completed successfully."


