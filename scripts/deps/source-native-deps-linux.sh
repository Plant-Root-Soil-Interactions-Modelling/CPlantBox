#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: scripts/deps/source-native-deps-linux.sh PREFIX

Build the native SuiteSparse/SUNDIALS dependencies needed by CPlantBox from
pinned source inputs and install headers/static libraries into PREFIX.

The output prefix is deleted and recreated. Do not pass a shared system prefix.
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

PREFIX=${1:-}
if [[ -z "${PREFIX}" ]]; then
  usage >&2
  exit 2
fi

PREFIX=$(python3 -c 'import os, sys; print(os.path.abspath(sys.argv[1]))' "${PREFIX}")
case "${PREFIX}" in
  /|/usr|/usr/local|/tmp|"${HOME:-__no_home__}")
    echo "Refusing to delete unsafe dependency prefix: ${PREFIX}" >&2
    exit 2
    ;;
esac

prefix_name=$(basename "${PREFIX}")
if [[ "${prefix_name}" != cplantbox-* ]]; then
  echo "Refusing to delete dependency prefix without cplantbox-* basename: ${PREFIX}" >&2
  echo "Use a dedicated prefix such as /tmp/cplantbox-native-deps." >&2
  exit 2
fi

JOBS=${CPB_DEPS_JOBS:-}
if [[ -z "${JOBS}" ]]; then
  if command -v nproc >/dev/null 2>&1; then
    JOBS=$(nproc)
  elif command -v sysctl >/dev/null 2>&1; then
    JOBS=$(sysctl -n hw.ncpu)
  else
    JOBS=2
  fi
fi

WORKDIR=${CPB_DEPS_WORKDIR:-}
if [[ -z "${WORKDIR}" ]]; then
  WORKDIR=$(mktemp -d "${TMPDIR:-/tmp}/cplantbox-native-deps.XXXXXXXXXX")
fi
WORKDIR=$(python3 -c 'import os, sys; print(os.path.abspath(sys.argv[1]))' "${WORKDIR}")
mkdir -p "${WORKDIR}"

SUNDIALS_VERSION=4.0.2
SUNDIALS_URL="https://github.com/LLNL/sundials/releases/download/v${SUNDIALS_VERSION}/sundials-${SUNDIALS_VERSION}.tar.gz"
SUNDIALS_SHA256="6656d6938aed9142e61a001b1ed9f4ee4f7eaf003613bf5a887e98a85904d375"

SUITESPARSE_VERSION=5.3.0
SUITESPARSE_URL="https://github.com/DrTimothyAldenDavis/SuiteSparse.git"
SUITESPARSE_TAG="v${SUITESPARSE_VERSION}"
SUITESPARSE_COMMIT="e927f7a3fc82339755482e553df37d932ff30083"

log() {
  printf '\n==> %s\n' "$*"
}

fetch() {
  local url=$1
  local output=$2
  if [[ -f "${output}" ]]; then
    return
  fi
  if command -v curl >/dev/null 2>&1; then
    curl -L --fail --retry 3 --retry-delay 2 -o "${output}" "${url}"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "${output}" "${url}"
  else
    echo "curl or wget is required to fetch ${url}" >&2
    exit 1
  fi
}

verify_sha256() {
  local expected=$1
  local file=$2
  python3 - "$expected" "$file" <<'PY'
import hashlib
import pathlib
import sys

expected = sys.argv[1]
path = pathlib.Path(sys.argv[2])
actual = hashlib.sha256(path.read_bytes()).hexdigest()
if actual != expected:
    raise SystemExit(f"SHA256 mismatch for {path}: expected {expected}, got {actual}")
print(f"{actual}  {path}")
PY
}

make_suite_sparse() {
  make -C "$1" static \
    AUTOCC=no \
    CC="${CC:-cc}" \
    CFLAGS="${CFLAGS:--O2 -fPIC}" \
    BLAS= \
    LAPACK= \
    CUDA=no
}

log "Preparing clean dependency prefix ${PREFIX}"
rm -rf "${PREFIX}"
mkdir -p "${PREFIX}/include" "${PREFIX}/lib" "${PREFIX}/share/licenses/suitesparse" "${PREFIX}/share/licenses/sundials"

log "Fetching SUNDIALS ${SUNDIALS_VERSION}"
SUNDIALS_TARBALL="${WORKDIR}/sundials-${SUNDIALS_VERSION}.tar.gz"
fetch "${SUNDIALS_URL}" "${SUNDIALS_TARBALL}"
verify_sha256 "${SUNDIALS_SHA256}" "${SUNDIALS_TARBALL}"
rm -rf "${WORKDIR}/sundials-${SUNDIALS_VERSION}"
tar -xzf "${SUNDIALS_TARBALL}" -C "${WORKDIR}"
SUNDIALS_SOURCE="${WORKDIR}/sundials-${SUNDIALS_VERSION}"

log "Fetching SuiteSparse ${SUITESPARSE_VERSION}"
SUITESPARSE_SOURCE="${WORKDIR}/SuiteSparse-${SUITESPARSE_VERSION}"
if [[ ! -d "${SUITESPARSE_SOURCE}/.git" ]]; then
  rm -rf "${SUITESPARSE_SOURCE}"
  git -c advice.detachedHead=false clone --depth 1 --branch "${SUITESPARSE_TAG}" "${SUITESPARSE_URL}" "${SUITESPARSE_SOURCE}"
fi
(
  cd "${SUITESPARSE_SOURCE}"
  git fetch --depth 1 origin "${SUITESPARSE_TAG}" >/dev/null 2>&1 || true
  git checkout --quiet "${SUITESPARSE_COMMIT}"
  actual_commit=$(git rev-parse HEAD)
  if [[ "${actual_commit}" != "${SUITESPARSE_COMMIT}" ]]; then
    echo "SuiteSparse commit mismatch: expected ${SUITESPARSE_COMMIT}, got ${actual_commit}" >&2
    exit 1
  fi
)

log "Building SuiteSparse static libraries"
for package in SuiteSparse_config AMD BTF COLAMD KLU; do
  make_suite_sparse "${SUITESPARSE_SOURCE}/${package}"
done

cp "${SUITESPARSE_SOURCE}/SuiteSparse_config/libsuitesparseconfig.a" "${PREFIX}/lib/"
cp "${SUITESPARSE_SOURCE}/AMD/Lib/libamd.a" "${PREFIX}/lib/"
cp "${SUITESPARSE_SOURCE}/BTF/Lib/libbtf.a" "${PREFIX}/lib/"
cp "${SUITESPARSE_SOURCE}/COLAMD/Lib/libcolamd.a" "${PREFIX}/lib/"
cp "${SUITESPARSE_SOURCE}/KLU/Lib/libklu.a" "${PREFIX}/lib/"

cp "${SUITESPARSE_SOURCE}/SuiteSparse_config/SuiteSparse_config.h" "${PREFIX}/include/"
cp "${SUITESPARSE_SOURCE}/AMD/Include/amd.h" "${PREFIX}/include/"
cp "${SUITESPARSE_SOURCE}/BTF/Include/btf.h" "${PREFIX}/include/"
cp "${SUITESPARSE_SOURCE}/COLAMD/Include/colamd.h" "${PREFIX}/include/"
cp "${SUITESPARSE_SOURCE}/KLU/Include/klu.h" "${PREFIX}/include/"

cp "${SUITESPARSE_SOURCE}/LICENSE.txt" "${PREFIX}/share/licenses/suitesparse/"
cp "${SUITESPARSE_SOURCE}/CONTRIBUTOR-LICENSE.txt" "${PREFIX}/share/licenses/suitesparse/"
for package in AMD BTF COLAMD KLU; do
  cp "${SUITESPARSE_SOURCE}/${package}/Doc/License.txt" "${PREFIX}/share/licenses/suitesparse/${package}-License.txt"
done

log "Patching SUNDIALS ${SUNDIALS_VERSION} for CMake 4 compatibility"
python3 - "${SUNDIALS_SOURCE}/config/SundialsKLU.cmake" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
text = path.read_text()
text = text.replace("CMAKE_MINIMUM_REQUIRED(VERSION 2.4)", "CMAKE_MINIMUM_REQUIRED(VERSION 3.5)")
path.write_text(text)
PY

log "Building SUNDIALS static libraries"
SUNDIALS_BUILD="${WORKDIR}/sundials-${SUNDIALS_VERSION}-build"
rm -rf "${SUNDIALS_BUILD}"
cmake -S "${SUNDIALS_SOURCE}" -B "${SUNDIALS_BUILD}" \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DCMAKE_C_FLAGS="${CFLAGS:--O2 -fPIC}" \
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_STATIC_LIBS=ON \
  -DBUILD_ARKODE=ON \
  -DBUILD_CVODE=ON \
  -DBUILD_CVODES=OFF \
  -DBUILD_IDA=OFF \
  -DBUILD_IDAS=OFF \
  -DBUILD_KINSOL=OFF \
  -DEXAMPLES_ENABLE_C=OFF \
  -DEXAMPLES_INSTALL=OFF \
  -DKLU_ENABLE=ON \
  -DKLU_INCLUDE_DIR="${PREFIX}/include" \
  -DKLU_LIBRARY="${PREFIX}/lib/libklu.a" \
  -DSUNDIALS_INDEX_SIZE=64
cmake --build "${SUNDIALS_BUILD}" -j"${JOBS}"
cmake --install "${SUNDIALS_BUILD}"

cp "${SUNDIALS_SOURCE}/LICENSE" "${PREFIX}/share/licenses/sundials/"
cp "${SUNDIALS_SOURCE}/NOTICE" "${PREFIX}/share/licenses/sundials/"

cat > "${PREFIX}/share/cplantbox-native-deps-provenance.txt" <<EOF
CPlantBox source-built native dependencies

SUNDIALS
  version: ${SUNDIALS_VERSION}
  source: ${SUNDIALS_URL}
  sha256: ${SUNDIALS_SHA256}
  build: static, KLU enabled, SUNDIALS_INDEX_SIZE=64

SuiteSparse
  version: ${SUITESPARSE_VERSION}
  source: ${SUITESPARSE_URL}
  tag: ${SUITESPARSE_TAG}
  commit: ${SUITESPARSE_COMMIT}
  components: SuiteSparse_config, AMD, BTF, COLAMD, KLU
  build: static, fPIC
EOF

log "Validating dependency prefix"
required_files=(
  include/klu.h
  include/amd.h
  include/btf.h
  include/colamd.h
  include/SuiteSparse_config.h
  include/sundials/sundials_config.h
  include/cvode/cvode.h
  include/nvector/nvector_serial.h
  include/sunlinsol/sunlinsol_klu.h
  include/sunmatrix/sunmatrix_sparse.h
  lib/libklu.a
  lib/libamd.a
  lib/libbtf.a
  lib/libcolamd.a
  lib/libsuitesparseconfig.a
  lib/libsundials_arkode.a
  lib/libsundials_cvode.a
  lib/libsundials_nvecserial.a
  lib/libsundials_sunlinsolklu.a
  lib/libsundials_sunmatrixsparse.a
  share/licenses/sundials/LICENSE
  share/licenses/sundials/NOTICE
  share/licenses/suitesparse/LICENSE.txt
  share/cplantbox-native-deps-provenance.txt
)
for file in "${required_files[@]}"; do
  if [[ ! -f "${PREFIX}/${file}" ]]; then
    echo "Missing expected dependency output: ${PREFIX}/${file}" >&2
    exit 1
  fi
done

log "Native dependency prefix ready: ${PREFIX}"
find "${PREFIX}/lib" -maxdepth 1 -type f -name '*.a' -print | sort
