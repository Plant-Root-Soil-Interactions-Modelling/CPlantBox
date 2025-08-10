#!/usr/bin/env bash
set -euo pipefail

# macOS (arm64/x86_64) developer helper: build editable install using system
# SUNDIALS and SuiteSparse from Homebrew, then run a headless smoke test.

if ! command -v brew >/dev/null 2>&1; then
  echo "Homebrew not found. Install from https://brew.sh and then:\n  brew install sundials suite-sparse" >&2
  exit 2
fi

BREW_PREFIX=$(brew --prefix)
need_pkgs=()
brew ls --versions sundials >/dev/null 2>&1 || need_pkgs+=(sundials)
brew ls --versions suite-sparse >/dev/null 2>&1 || need_pkgs+=(suite-sparse)
if ((${#need_pkgs[@]})); then
  echo "Installing required packages: ${need_pkgs[*]}"
  brew install "${need_pkgs[@]}"
fi

# Ensure Python venv
if [[ ! -d .venv ]]; then
  python3 -m venv .venv
fi
source .venv/bin/activate
python -m pip -q install -U pip

# Prefer system libs; disable bundled for mac dev builds
export CMAKE_PREFIX_PATH="${BREW_PREFIX}:${CMAKE_PREFIX_PATH:-}"
export CMAKE_ARGS="-DUSE_SYSTEM_SUITESPARSE=ON -DBUNDLE_SUITESPARSE=OFF -DUSE_SYSTEM_SUNDIALS=ON -DBUNDLE_SUNDIALS=OFF ${CMAKE_ARGS:-}"

python -m pip -q install -e .

python - <<'PY'
import os
import plantbox as pb
print("Imported plantbox", getattr(pb, "__version__", "?"))
rs = pb.Plant()
try:
    rs.readParameters(os.path.join("modelparameter","structural","rootsystem","Anagallis_femina_Leitner_2010.xml"), verbose=False)
except Exception:
    rs.readParameters(os.path.join(pb.data_path(), "structural","rootsystem","Anagallis_femina_Leitner_2010.xml"), verbose=False)
rs.initialize(False)
rs.simulate(1.0, False)
print("Editable install import/simulate OK (macOS)")
PY


