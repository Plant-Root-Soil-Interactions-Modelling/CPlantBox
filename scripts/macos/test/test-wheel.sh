#!/usr/bin/env bash
set -euo pipefail

# Usage: scripts/macos/test/test-wheel.sh [WHEEL_PATH]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$REPO_ROOT"

WHEEL="${1:-}"
if [[ -z "$WHEEL" ]]; then
  ARCH=$(uname -m)
  case "$ARCH" in
    arm64) DIR="wheelhouse/macos"; SELECTOR="*arm64.whl" ;;
    x86_64) DIR="wheelhouse/macos"; SELECTOR="*x86_64.whl" ;;
    *) echo "Unsupported macOS arch: $ARCH" >&2; exit 2 ;;
  esac
  WHEEL=$(ls -1t "$DIR"/$SELECTOR | head -n1)
fi
[[ -f "$WHEEL" ]] || { echo "Wheel not found: $WHEEL" >&2; exit 2; }

python3 -m venv /tmp/venv_cpb_mac --clear
source /tmp/venv_cpb_mac/bin/activate
python3 -m pip -q install -U pip
# Prefer a stable vtk version; allow override via VTK_VERSION env
VTK_VERSION=${VTK_VERSION:-9.2.6}
python3 -m pip -q install "vtk==${VTK_VERSION}" || python3 -m pip -q install vtk
python3 -m pip -q install "$WHEEL"

export PYTHONUNBUFFERED=1
export LIBGL_ALWAYS_SOFTWARE=1
export PYTHONPATH="$REPO_ROOT/test"
export CPB_GOLDEN_OS=macos

python3 scripts/common/golden_test.py

#!/usr/bin/env bash
set -euo pipefail

# Test a macOS wheel by installing it into a fresh venv and running:
# 1) minimal import/simulate smoke
# 2) headless golden image comparison
# Usage: scripts/macos/test/test-wheel.sh [WHEEL_PATH]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$REPO_ROOT"

source "$REPO_ROOT/scripts/common/env.sh"
source "$REPO_ROOT/scripts/common/logging.sh"

WHEEL_PATH_INPUT=${1:-}

if [[ -z "$WHEEL_PATH_INPUT" ]]; then
  WHEEL_DIR="wheelhouse/macos"
  if [[ ! -d "$WHEEL_DIR" ]]; then
    die "No wheelhouse dir found at $WHEEL_DIR. Run macOS build first."
  fi
  WHEEL_PATH_INPUT=$(ls -1t "$WHEEL_DIR"/*.whl | head -n1)
fi

if [[ ! -f "$WHEEL_PATH_INPUT" ]]; then
  die "Wheel not found: $WHEEL_PATH_INPUT"
fi

# Convert to absolute path for clarity
WHEEL_ABS=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$WHEEL_PATH_INPUT")

log_info "[macOS-test] Using wheel: $WHEEL_ABS"

# Reuse the same inner logic as a dedicated runner script
"$REPO_ROOT/scripts/macos/test/run_golden.sh" "$WHEEL_ABS"

log_info "[macOS-test] Completed."


