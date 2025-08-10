#!/usr/bin/env bash
set -euo pipefail

# Usage: run_golden_in_container.sh [x86_64|aarch64] <wheel_path>

ARCH=${1:-x86_64}
WHEEL=${2:-}
[[ -n "$WHEEL" && -f "$WHEEL" ]] || { echo "Usage: $0 [x86_64|aarch64] <wheel_path>" >&2; exit 2; }

case "$ARCH" in
  x86_64) IMG=cplantbox-manylinux-test-x86_64; PY=/opt/python/cp310-cp310/bin/python;;
  aarch64) IMG=cplantbox-manylinux-test-aarch64; PY=/opt/python/cp310-cp310/bin/python;;
  *) echo "Unsupported arch: $ARCH" >&2; exit 2;;
esac

docker run --rm --platform linux/${ARCH} -v "$PWD":/project -w /project "$IMG" \
  bash -lc "set -euo pipefail; ${PY} -m pip install -U pip pytest vtk==9.2.6; ${PY} -m pip install '$WHEEL'; export PYTHONUNBUFFERED=1; export LIBGL_ALWAYS_SOFTWARE=1; export MESA_LOADER_DRIVER_OVERRIDE=llvmpipe; export PYTHONPATH=/project/test; export CPB_GOLDEN_OS=linux; Xvfb :99 -screen 0 1280x1024x24 +extension RANDR >/tmp/xvfb.log 2>&1 & sleep 1; export DISPLAY=:99; ${PY} scripts/common/golden_test.py"


