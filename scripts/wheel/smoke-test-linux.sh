#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: scripts/wheel/smoke-test-linux.sh

Build a CPlantBox wheel with source-built native dependencies, install it into
an isolated virtual environment, and run scripts/wheel/smoke_test.py.

Environment overrides:
  CPB_NATIVE_DEPS_PREFIX   Native dependency prefix [default: /tmp/cplantbox-native-deps]
  CPB_DEPS_WORKDIR         Source dependency working directory
  CPB_DEPS_JOBS            Source dependency parallelism
  CPB_WHEELHOUSE           Wheel output directory [default: /tmp/cplantbox-wheelhouse]
  CPB_WHEEL_BUILD_VENV     Build virtualenv [default: /tmp/cplantbox-wheel-build]
  CPB_WHEEL_TEST_VENV      Test virtualenv [default: /tmp/cplantbox-wheel-test]
  CPB_SKBUILD_DIR          scikit-build-core build dir [default: /tmp/cplantbox-skbuild]
  CMAKE_ARGS               Extra CMake arguments appended to the source-provider arguments

The helper deletes/recreates CPB_WHEELHOUSE, CPB_WHEEL_BUILD_VENV,
CPB_WHEEL_TEST_VENV, and CPB_SKBUILD_DIR. Use dedicated cplantbox-* paths.
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "${SCRIPT_DIR}/../.." && pwd)

DEPS_PREFIX=${CPB_NATIVE_DEPS_PREFIX:-/tmp/cplantbox-native-deps}
WHEELHOUSE=${CPB_WHEELHOUSE:-/tmp/cplantbox-wheelhouse}
BUILD_VENV=${CPB_WHEEL_BUILD_VENV:-/tmp/cplantbox-wheel-build}
TEST_VENV=${CPB_WHEEL_TEST_VENV:-/tmp/cplantbox-wheel-test}
SKBUILD_DIR=${CPB_SKBUILD_DIR:-/tmp/cplantbox-skbuild}

abs_path() {
  python3 -c 'import os, sys; print(os.path.abspath(sys.argv[1]))' "$1"
}

require_safe_output_path() {
  local path
  path=$(abs_path "$1")
  case "${path}" in
    /|/usr|/usr/local|/tmp|"${HOME:-__no_home__}")
      echo "Refusing to delete unsafe wheel smoke-test path: ${path}" >&2
      exit 2
      ;;
  esac
  if [[ "$(basename "${path}")" != cplantbox-* ]]; then
    echo "Refusing to delete path without cplantbox-* basename: ${path}" >&2
    exit 2
  fi
}

for path in "${WHEELHOUSE}" "${BUILD_VENV}" "${TEST_VENV}" "${SKBUILD_DIR}"; do
  require_safe_output_path "${path}"
done

WHEELHOUSE=$(abs_path "${WHEELHOUSE}")
BUILD_VENV=$(abs_path "${BUILD_VENV}")
TEST_VENV=$(abs_path "${TEST_VENV}")
SKBUILD_DIR=$(abs_path "${SKBUILD_DIR}")

log() {
  printf '\n==> %s\n' "$*"
}

log "Building source native dependencies"
"${REPO_ROOT}/scripts/deps/source-native-deps-linux.sh" "${DEPS_PREFIX}"

log "Creating wheel build environment"
rm -rf "${WHEELHOUSE}" "${BUILD_VENV}" "${TEST_VENV}" "${SKBUILD_DIR}"
mkdir -p "${WHEELHOUSE}"
python3 -m venv "${BUILD_VENV}"
# shellcheck disable=SC1091
source "${BUILD_VENV}/bin/activate"
python -m pip install --upgrade pip build

log "Building wheel into ${WHEELHOUSE}"
export CMAKE_ARGS="-DCPB_SUITESPARSE_PROVIDER=source -DCPB_SUNDIALS_PROVIDER=source -DCPB_NATIVE_DEPS_PREFIX=${DEPS_PREFIX} ${CMAKE_ARGS:-}"
python -m build --wheel --outdir "${WHEELHOUSE}" --config-setting="build-dir=${SKBUILD_DIR}" "${REPO_ROOT}"
deactivate

wheel_count=$(find "${WHEELHOUSE}" -maxdepth 1 -type f -name '*.whl' | wc -l | tr -d ' ')
if [[ "${wheel_count}" != "1" ]]; then
  echo "Expected exactly one wheel in ${WHEELHOUSE}, found ${wheel_count}" >&2
  find "${WHEELHOUSE}" -maxdepth 1 -type f -print >&2
  exit 1
fi
wheel=$(find "${WHEELHOUSE}" -maxdepth 1 -type f -name '*.whl' | sort | head -n 1)

log "Installing wheel into isolated test environment"
python3 -m venv "${TEST_VENV}"
# shellcheck disable=SC1091
source "${TEST_VENV}/bin/activate"
python -m pip install --upgrade pip
python -m pip install --no-deps "${wheel}"

log "Running installed-wheel smoke test"
cd /tmp
python "${REPO_ROOT}/scripts/wheel/smoke_test.py"
