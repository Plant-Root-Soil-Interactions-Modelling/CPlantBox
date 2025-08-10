#!/usr/bin/env bash
set -euo pipefail

# Build macOS wheels for CPython 3.9–3.12 using system SUNDIALS/SuiteSparse (Homebrew),
# then repair with delocate and run a headless smoke test.

# Ensure we always run from repo root and can restore it after temp cd's
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

if ! command -v brew >/dev/null 2>&1; then
  echo "Homebrew not found. Install from https://brew.sh" >&2
  exit 2
fi

BREW_PREFIX=$(brew --prefix)
python_bins=(
  $(ls -1 /usr/local/bin/python3.{9,10,11,12} 2>/dev/null || true)
  $(ls -1 /opt/homebrew/bin/python3.{9,10,11,12} 2>/dev/null || true)
)

if [[ ${#python_bins[@]} -eq 0 ]]; then
  echo "No CPython 3.9–3.12 found in standard locations. Install via pyenv or Homebrew." >&2
  exit 2
fi

pipx -q >/dev/null 2>&1 || true

need_pkgs=()
brew ls --versions sundials >/dev/null 2>&1 || need_pkgs+=(sundials)
brew ls --versions suite-sparse >/dev/null 2>&1 || need_pkgs+=(suite-sparse)
if ((${#need_pkgs[@]})); then
  echo "Installing required packages: ${need_pkgs[*]}"
  brew install "${need_pkgs[@]}"
fi

python3 -m pip install -U delocate build >/dev/null

OUT_WHL_DIR="wheelhouse/macos"
rm -rf "$OUT_WHL_DIR"
mkdir -p "$OUT_WHL_DIR"
export MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-11.0}"
for PY in "${python_bins[@]}"; do
  echo "[macOS] Building with $PY"
  # Build inside a temp venv to avoid PEP 668 on Homebrew Pythons
  TMPBUILD=$(mktemp -d)
  "$PY" -m venv "$TMPBUILD/venv"
  "$TMPBUILD/venv/bin/python" -m pip install -U pip build >/dev/null
  rm -rf dist/ _skbuild/
  export CMAKE_PREFIX_PATH="${BREW_PREFIX}:${CMAKE_PREFIX_PATH:-}"
  export CMAKE_ARGS="-DUSE_SYSTEM_SUITESPARSE=ON -DBUNDLE_SUITESPARSE=OFF -DUSE_SYSTEM_SUNDIALS=ON -DBUNDLE_SUNDIALS=OFF -DCMAKE_OSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET} ${CMAKE_ARGS:-}"
  "$TMPBUILD/venv/bin/python" -m build --wheel
  ls -lh dist/*.whl
  WHL=$(ls -1 dist/*.whl | head -n1)
  delocate-listdeps "$WHL" || true
  delocate-wheel -v -w "$OUT_WHL_DIR" "$WHL"

  echo "[macOS] Testing repaired wheel for $PY"
  TMPVENV=$(mktemp -d)
  "$PY" -m venv "$TMPVENV/venv"
  "$TMPVENV/venv/bin/python" -m pip install -U pip >/dev/null
  # Install from the repaired wheel directory; pip will select the correct wheel for this interpreter
  "$TMPVENV/venv/bin/python" -m pip install --no-index --find-links="$OUT_WHL_DIR" "cplantbox==2.1.0" >/dev/null
  cd "$(mktemp -d)"
  "$TMPVENV/venv/bin/python" - <<'PY'
import os
import plantbox as pb
rs = pb.Plant()
# Always resolve from packaged data to avoid silent non-throwing failures
p = os.path.join(pb.data_path(), "structural","rootsystem","Anagallis_femina_Leitner_2010.xml")
rs.readParameters(p, verbose=False)
rs.initialize(False)
rs.simulate(1.0, False)
print("macOS wheel import/simulate OK")
PY
  cd "$REPO_ROOT"
  rm -rf "$TMPVENV"
  rm -rf "$TMPBUILD"
done

echo "[macOS] Done. Artifacts in wheelhouse/"


