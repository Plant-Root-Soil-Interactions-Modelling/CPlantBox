#!/usr/bin/env bash
set -euo pipefail

# Install a provided wheel into isolated venvs and run both a minimal
# import/simulate smoke and a headless golden-image comparison.
# Usage: scripts/macos/test/run_golden.sh <wheel_path>

WHEEL_PATH=${1:-}
if [[ -z "$WHEEL_PATH" ]]; then
  echo "Usage: $0 <wheel_path>" >&2
  exit 2
fi

if [[ ! -f "$WHEEL_PATH" ]]; then
  echo "Wheel not found: $WHEEL_PATH" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$REPO_ROOT"

source "$REPO_ROOT/scripts/common/env.sh"
source "$REPO_ROOT/scripts/common/logging.sh"

log_info "[macOS-run] Creating venv for wheel smoke..."
SMOKE_DIR=$(mktemp -d 2>/dev/null || mktemp -d -t cplantbox-smoke)
python3 -m venv "$SMOKE_DIR/venv"
source "$SMOKE_DIR/venv/bin/activate"
python3 -m pip -q install -U pip
python3 -m pip -q install "$WHEEL_PATH"

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

TMPDIR1=$(mktemp -d 2>/dev/null || mktemp -d -t cplantbox-smoke)
pushd "$TMPDIR1" >/dev/null
export REPO_ROOT_ESCAPED="$REPO_ROOT"
python3 - <<'PY'
import os, sys, sysconfig
REPO_ROOT = os.environ.get("REPO_ROOT_ESCAPED", "")
site = sysconfig.get_paths()["platlib"]
# Ensure installed site takes precedence but keep stdlib
if REPO_ROOT:
    sys.path = [p for p in sys.path if not p.startswith(REPO_ROOT)]
if site not in sys.path:
    sys.path.insert(0, site)
import plantbox as pb
print("Imported plantbox", getattr(pb, "__version__", "?"))
print("Module file:", getattr(pb, "__file__", "?"))
rootsys = os.path.join(pb.data_path(), "structural", "rootsystem", "Anagallis_femina_Leitner_2010.xml")
print("Param exists:", os.path.exists(rootsys))
rs = pb.Plant()
rs.readParameters(rootsys, verbose=False)
rs.initialize(False)
rs.simulate(2.0, False)
print("Wheel import/simulate OK (macOS)")
PY
popd >/dev/null

deactivate
rm -rf "$SMOKE_DIR"

log_info "[macOS-run] Creating venv for golden test..."
TEST_DIR=$(mktemp -d 2>/dev/null || mktemp -d -t cplantbox-test)
python3 -m venv "$TEST_DIR/venv"
source "$TEST_DIR/venv/bin/activate"
python3 -m pip -q install -U pip pytest vtk
python3 -m pip -q install "$WHEEL_PATH"

# Use the macOS golden test runner that isolates sys.path and imports the wheel
python3 scripts/macos/test/golden_test.py

deactivate
rm -rf "$TEST_DIR"

log_info "[macOS-run] Golden smoke finished."


