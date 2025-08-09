#!/usr/bin/env bash
set -euo pipefail

# Run curated tests and the headless golden-image test on macOS host.
# Assumes an active Python environment where plantbox is installed (editable or wheel)
# and optional VTK is available for offscreen rendering.

if [ -d .venv ]; then
  source .venv/bin/activate
fi

echo "Python: $(python3 -V || true)"
echo "Which python: $(which python3 || true)"

# Prefer deterministic behavior in CI/local runs
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

echo "Generating golden image via tut_1_3_headless.py..."
python3 tut_1_3_headless.py

echo "Attempting wheel build and wheel-based smoke (similar to Ubuntu script)..."
if python3 -m pip -q install build && rm -rf dist/ _skbuild/ && SETUPTOOLS_SCM_PRETEND_VERSION="$(git describe --tags --abbrev=0 2>/dev/null || echo 2.1.0)" python3 -m build --wheel >/dev/null 2>&1; then
  echo "Wheel built; testing installed package..."
  python3 -m venv /tmp/venv && . /tmp/venv/bin/activate
  python3 -m pip -q install --upgrade pip
  WHEEL_PATH=$(ls -1t dist/*.whl | head -n1)
  python3 -m pip -q install "$WHEEL_PATH"
  python3 - <<PY
import plantbox as pb
rs = pb.Plant()
rs.readParameters("modelparameter/structural/rootsystem/Anagallis_femina_Leitner_2010.xml", verbose=False)
rs.initialize(False)
rs.simulate(5.0, False)
print("Wheel import/simulate OK (macOS)")
PY
  deactivate
else
  echo "Wheel build failed; falling back to source-tree curated tests..."
  export PYTHONPATH="$(pwd):$(pwd)/src:${PYTHONPATH:-}"
  pushd test >/dev/null
  pytest -q \
    test_organ.py \
    test_organparameter.py \
    test_rootparameter.py \
    test_leafparameter.py \
    test_leaf.py \
    test_hydraulic_parameters.py
  popd >/dev/null
fi

echo "Running golden image comparison test..."
python3 -m pytest -q test/test_golden_headless.py

echo "All macOS tests completed."


