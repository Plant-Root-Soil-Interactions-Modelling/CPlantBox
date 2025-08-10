#!/usr/bin/env bash
set -euo pipefail
set -x

WHEEL_PATH_REL=${1:-}
if [[ -z "$WHEEL_PATH_REL" ]]; then
  echo "Usage: $0 <wheel_path_relative_to_repo_root>" >&2
  exit 2
fi

if [[ ! -f "$WHEEL_PATH_REL" ]]; then
  echo "Wheel not found: $WHEEL_PATH_REL" >&2
  exit 2
fi

python3 -m venv /tmp/venv --system-site-packages
. /tmp/venv/bin/activate
python3 -m pip -q install -U pip
python3 -m pip -q install "$WHEEL_PATH_REL"

export PYTHONUNBUFFERED=1
export LIBGL_ALWAYS_SOFTWARE=1
export MESA_LOADER_DRIVER_OVERRIDE=llvmpipe
export PYTHONPATH=/src/test
export CPB_GOLDEN_OS=linux

timeout 180s xvfb-run -a -s '-screen 0 1280x1024x24' \
  python3 scripts/common/golden_test.py


