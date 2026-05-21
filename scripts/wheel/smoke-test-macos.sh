#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: scripts/wheel/smoke-test-macos.sh

Build a macOS CPlantBox wheel with source-built native dependencies, install it
into an isolated virtual environment, audit it with delocate-listdeps/otool, and
run scripts/wheel/smoke_test.py outside the source tree.

Environment overrides:
  CPB_NATIVE_DEPS_PREFIX   Native dependency prefix [default: /tmp/cplantbox-native-deps-macos]
  CPB_DEPS_WORKDIR         Source dependency working directory
  CPB_DEPS_JOBS            Source dependency parallelism
  CPB_WHEELHOUSE           Wheel output directory [default: /tmp/cplantbox-wheelhouse-macos]
  CPB_WHEEL_BUILD_VENV     Build virtualenv [default: /tmp/cplantbox-wheel-build-macos]
  CPB_WHEEL_TEST_VENV      Test virtualenv [default: /tmp/cplantbox-wheel-test-macos]
  CPB_SKBUILD_DIR          scikit-build-core build dir [default: /tmp/cplantbox-skbuild-macos]
  CPB_PYTHON               Python interpreter for build/test venvs [default: python3]
  CPB_BUILD_NATIVE_DEPS    Build native deps before wheel build: 1|0 [default: 1]
  CMAKE_ARGS               Extra CMake arguments appended to the source-provider arguments

The helper deletes/recreates CPB_WHEELHOUSE, CPB_WHEEL_BUILD_VENV,
CPB_WHEEL_TEST_VENV, and CPB_SKBUILD_DIR. Use dedicated cplantbox-* paths.
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "smoke-test-macos.sh must run on macOS/Darwin" >&2
  exit 2
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "${SCRIPT_DIR}/../.." && pwd)

DEPS_PREFIX=${CPB_NATIVE_DEPS_PREFIX:-/tmp/cplantbox-native-deps-macos}
WHEELHOUSE=${CPB_WHEELHOUSE:-/tmp/cplantbox-wheelhouse-macos}
BUILD_VENV=${CPB_WHEEL_BUILD_VENV:-/tmp/cplantbox-wheel-build-macos}
TEST_VENV=${CPB_WHEEL_TEST_VENV:-/tmp/cplantbox-wheel-test-macos}
SKBUILD_DIR=${CPB_SKBUILD_DIR:-/tmp/cplantbox-skbuild-macos}
PYTHON=${CPB_PYTHON:-python3}
BUILD_NATIVE_DEPS=${CPB_BUILD_NATIVE_DEPS:-1}

abs_path() {
  python3 -c 'import os, sys; print(os.path.abspath(sys.argv[1]))' "$1"
}

require_safe_output_path() {
  local path
  path=$(abs_path "$1")
  case "${path}" in
    /|/usr|/usr/local|/tmp|/private/tmp|"${HOME:-__no_home__}")
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

DEPS_PREFIX=$(abs_path "${DEPS_PREFIX}")
WHEELHOUSE=$(abs_path "${WHEELHOUSE}")
BUILD_VENV=$(abs_path "${BUILD_VENV}")
TEST_VENV=$(abs_path "${TEST_VENV}")
SKBUILD_DIR=$(abs_path "${SKBUILD_DIR}")

log() {
  printf '\n==> %s\n' "$*"
}

export MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-11.0}
export CFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} ${CFLAGS:-}"
export CXXFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} ${CXXFLAGS:-}"

case "${BUILD_NATIVE_DEPS}" in
  0|1) ;;
  *) echo "CPB_BUILD_NATIVE_DEPS must be 0 or 1, got: ${BUILD_NATIVE_DEPS}" >&2; exit 2 ;;
esac

if [[ "${BUILD_NATIVE_DEPS}" == "1" ]]; then
  log "Building source native dependencies"
  "${REPO_ROOT}/scripts/deps/source-native-deps-macos.sh" "${DEPS_PREFIX}"
elif [[ ! -f "${DEPS_PREFIX}/share/cplantbox-native-deps-provenance.txt" ]]; then
  echo "CPB_BUILD_NATIVE_DEPS=0 but native dependency prefix is missing provenance: ${DEPS_PREFIX}" >&2
  exit 2
fi

log "Creating wheel build environment with ${PYTHON}"
rm -rf "${WHEELHOUSE}" "${BUILD_VENV}" "${TEST_VENV}" "${SKBUILD_DIR}"
mkdir -p "${WHEELHOUSE}"
"${PYTHON}" -m venv "${BUILD_VENV}"
# shellcheck disable=SC1091
source "${BUILD_VENV}/bin/activate"
python -m pip install --upgrade pip build delocate

log "Building wheel into ${WHEELHOUSE}"
export CMAKE_ARGS="-DCPB_SUITESPARSE_PROVIDER=source -DCPB_SUNDIALS_PROVIDER=source -DCPB_NATIVE_DEPS_PREFIX=${DEPS_PREFIX} -DCMAKE_OSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET} ${CMAKE_ARGS:-}"
python -m build --wheel --outdir "${WHEELHOUSE}" --config-setting="build-dir=${SKBUILD_DIR}" "${REPO_ROOT}"

wheel_count=$(find "${WHEELHOUSE}" -maxdepth 1 -type f -name '*.whl' | wc -l | tr -d ' ')
if [[ "${wheel_count}" != "1" ]]; then
  echo "Expected exactly one wheel in ${WHEELHOUSE}, found ${wheel_count}" >&2
  find "${WHEELHOUSE}" -maxdepth 1 -type f -print >&2
  exit 1
fi
wheel=$(find "${WHEELHOUSE}" -maxdepth 1 -type f -name '*.whl' | sort | head -n 1)

log "Auditing macOS wheel dependencies"
DELOCATE_LISTDEPS="${BUILD_VENV}/bin/delocate-listdeps" \
  CPB_EXPECTED_WHEEL_COUNT=1 \
  "${REPO_ROOT}/scripts/ci/audit-wheelhouse-macos.sh" "${WHEELHOUSE}"

deactivate

log "Installing wheel into isolated test environment"
"${PYTHON}" -m venv "${TEST_VENV}"
# shellcheck disable=SC1091
source "${TEST_VENV}/bin/activate"
python -m pip install --upgrade pip
python -m pip install --no-deps "${wheel}"

log "Running installed-wheel smoke test"
cd /tmp
python "${REPO_ROOT}/scripts/wheel/smoke_test.py"
