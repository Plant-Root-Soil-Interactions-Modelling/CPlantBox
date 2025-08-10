#!/usr/bin/env bash
set -euo pipefail

# Test a manylinux-repaired wheel inside a matching manylinux container, headless.
# Usage: scripts/wheels/linux/manylinux/test/test-wheel.sh [x86_64|aarch64] [WHEEL_PATH]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../../../" && pwd)"
cd "$REPO_ROOT"

ARCH="${1:-x86_64}"
case "$ARCH" in
  x86_64|aarch64) ;;
  *) echo "Unsupported arch: $ARCH" >&2; exit 2;;
esac

WHEEL_ARG="${2:-}"
if [[ -z "$WHEEL_ARG" ]]; then
  DIR="wheelhouse/linux/manylinux2014_${ARCH}/"
  WHEEL_ARG=$(ls -1t "$DIR"/*.whl | head -n1)
fi
[[ -f "$WHEEL_ARG" ]] || { echo "Wheel not found: $WHEEL_ARG" >&2; exit 2; }

if [[ "$ARCH" == "x86_64" ]]; then
  IMG="quay.io/pypa/manylinux2014_x86_64"
  DF="$SCRIPT_DIR/Dockerfile_x86_64"
  PY="/opt/python/cp310-cp310/bin/python"
else
  IMG="quay.io/pypa/manylinux2014_aarch64"
  DF="$SCRIPT_DIR/Dockerfile_aarch64"
  PY="/opt/python/cp310-cp310/bin/python"
fi

docker buildx build --platform linux/${ARCH} -f "$DF" -t cplantbox-manylinux-test-${ARCH} --load .

# Run golden in manylinux container; mount repo at /project
docker run --rm --platform linux/${ARCH} -v "$REPO_ROOT":/project -w /project cplantbox-manylinux-test-${ARCH} \
  bash -lc "set -euo pipefail; ${PY} -m pip install '$WHEEL_ARG'; export PYTHONPATH=/project/test; xvfb-run -a -s '-screen 0 1280x1024x24' ${PY} scripts/wheels/linux/golden_smoke.py"


