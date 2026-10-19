#!/usr/bin/env bash
set -euo pipefail

# Windows (MSVC) counterpart of scripts/deps/source-native-deps.sh.
#
# Builds the same pinned SuiteSparse/SUNDIALS dependencies as static MSVC
# libraries and installs headers/import libraries into PREFIX. Runs under the
# Git Bash that GitHub windows runners provide; the actual compilation is done
# by CMake with the Visual Studio generator, so no MSVC environment setup is
# needed here.
#
# Differences from the POSIX script, all forced by MSVC:
#   * SuiteSparse 5.x only ships Makefiles, so the packages are built through
#     scripts/deps/suitesparse-msvc/CMakeLists.txt, which reproduces the
#     Makefile variant matrix exactly.
#   * Static libraries are named klu.lib etc. rather than libklu.a; the
#     find_library calls in cmake/CPlantBoxPiafMunch.cmake accept both.

usage() {
	cat <<'USAGE'
Usage: scripts/deps/source-native-deps-windows.sh PREFIX

Build the native SuiteSparse/SUNDIALS dependencies needed by CPlantBox from
pinned source inputs with MSVC and install headers/static libraries into
PREFIX.

The output prefix is deleted and recreated. Do not pass a shared system prefix.
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
	usage
	exit 0
fi

PYTHON=python3
command -v python3 >/dev/null 2>&1 || PYTHON=python

PREFIX=${1:-}
if [[ -z "${PREFIX}" ]]; then
	usage >&2
	exit 2
fi

PREFIX=$("${PYTHON}" -c 'import os, sys; print(os.path.abspath(sys.argv[1]).replace("\\", "/"))' "${PREFIX}")
prefix_name=$(basename "${PREFIX}")
if [[ "${prefix_name}" != cplantbox-* ]]; then
	echo "Refusing to delete dependency prefix without cplantbox-* basename: ${PREFIX}" >&2
	echo "Use a dedicated prefix such as C:/cplantbox-native-deps." >&2
	exit 2
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

log() {
	echo "[native-deps-windows] $*"
}

verify_sha256() {
	local expected=$1
	local path=$2
	"${PYTHON}" - "$expected" "$path" <<'PY'
import hashlib
import pathlib
import sys

expected = sys.argv[1]
path = pathlib.Path(sys.argv[2])
actual = hashlib.sha256(path.read_bytes()).hexdigest()
if actual != expected:
    raise SystemExit(f"sha256 mismatch for {path}: expected {expected}, got {actual}")
PY
}

WORKDIR=${CPB_DEPS_WORKDIR:-}
if [[ -z "${WORKDIR}" ]]; then
	WORKDIR=$(mktemp -d "${TMPDIR:-/tmp}/cplantbox-native-deps.XXXXXXXXXX")
fi
WORKDIR=$("${PYTHON}" -c 'import os, sys; print(os.path.abspath(sys.argv[1]).replace("\\", "/"))' "${WORKDIR}")
mkdir -p "${WORKDIR}"

JOBS=${CPB_DEPS_JOBS:-4}

SUNDIALS_VERSION=4.0.2
SUNDIALS_URL="https://github.com/LLNL/sundials/releases/download/v${SUNDIALS_VERSION}/sundials-${SUNDIALS_VERSION}.tar.gz"
SUNDIALS_SHA256="6656d6938aed9142e61a001b1ed9f4ee4f7eaf003613bf5a887e98a85904d375"

SUITESPARSE_VERSION=5.3.0
SUITESPARSE_URL="https://github.com/DrTimothyAldenDavis/SuiteSparse.git"
SUITESPARSE_TAG="v${SUITESPARSE_VERSION}"
SUITESPARSE_COMMIT="e927f7a3fc82339755482e553df37d932ff30083"

CMAKE_GENERATOR_ARGS=(-A x64)

log "Preparing clean dependency prefix ${PREFIX}"
rm -rf "${PREFIX}"
mkdir -p "${PREFIX}/include" "${PREFIX}/lib" "${PREFIX}/share/licenses/suitesparse" "${PREFIX}/share/licenses/sundials"

log "Fetching SUNDIALS ${SUNDIALS_VERSION}"
SUNDIALS_TARBALL="${WORKDIR}/sundials-${SUNDIALS_VERSION}.tar.gz"
if [[ ! -f "${SUNDIALS_TARBALL}" ]]; then
	curl -fsSL "${SUNDIALS_URL}" -o "${SUNDIALS_TARBALL}"
fi
verify_sha256 "${SUNDIALS_SHA256}" "${SUNDIALS_TARBALL}"
rm -rf "${WORKDIR}/sundials-${SUNDIALS_VERSION}"
# --force-local: GNU tar would otherwise treat the drive-letter colon in the
# Windows path as a remote-host separator
tar --force-local -xzf "${SUNDIALS_TARBALL}" -C "${WORKDIR}"
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

log "Building SuiteSparse static libraries with MSVC"
SUITESPARSE_BUILD="${WORKDIR}/suitesparse-msvc-build"
rm -rf "${SUITESPARSE_BUILD}"
cmake -S "${SCRIPT_DIR}/suitesparse-msvc" -B "${SUITESPARSE_BUILD}" \
	"${CMAKE_GENERATOR_ARGS[@]}" \
	-DSUITESPARSE_SOURCE_DIR="${SUITESPARSE_SOURCE}"
cmake --build "${SUITESPARSE_BUILD}" --config Release -j "${JOBS}"

for lib in suitesparseconfig amd btf colamd klu; do
	cp "${SUITESPARSE_BUILD}/Release/${lib}.lib" "${PREFIX}/lib/"
done
# SUNDIALS' FindKLU prepends "lib" to CMAKE_FIND_LIBRARY_PREFIXES in a way that
# loses the empty prefix, so on Windows it can only find lib-prefixed names
# (except suitesparseconfig, which it special-cases without a prefix). Install
# lib-prefixed copies alongside the plain ones so both consumers are satisfied.
for lib in amd btf colamd klu; do
	cp "${SUITESPARSE_BUILD}/Release/${lib}.lib" "${PREFIX}/lib/lib${lib}.lib"
done

cp "${SUITESPARSE_SOURCE}/SuiteSparse_config/SuiteSparse_config.h" "${PREFIX}/include/"
cp "${SUITESPARSE_SOURCE}/AMD/Include/amd.h" "${PREFIX}/include/"
cp "${SUITESPARSE_SOURCE}/BTF/Include/btf.h" "${PREFIX}/include/"
cp "${SUITESPARSE_SOURCE}/COLAMD/Include/colamd.h" "${PREFIX}/include/"
cp "${SUITESPARSE_SOURCE}/KLU/Include/klu.h" "${PREFIX}/include/"

log "DIAG: prefix contents"
ls -l "${PREFIX}/lib" "${PREFIX}/include" || true

log "DIAG: find_library probe"
PROBE_DIR="${WORKDIR}/find-probe"
rm -rf "${PROBE_DIR}"
mkdir -p "${PROBE_DIR}"
cat >"${PROBE_DIR}/CMakeLists.txt" <<EOF
cmake_minimum_required(VERSION 3.18)
project(probe C)
set(d "${PREFIX}/lib")
find_library(P_SHORT amd \${d} NO_DEFAULT_PATH)
find_library(P_KEYWORD NAMES amd libamd PATHS \${d} NO_DEFAULT_PATH)
message(STATUS "PROBE short=\${P_SHORT}")
message(STATUS "PROBE keyword=\${P_KEYWORD}")
message(STATUS "PROBE suffixes=\${CMAKE_FIND_LIBRARY_SUFFIXES} prefixes=\${CMAKE_FIND_LIBRARY_PREFIXES}")
EOF
cmake -S "${PROBE_DIR}" -B "${PROBE_DIR}/build" "${CMAKE_GENERATOR_ARGS[@]}" || true

cp "${SUITESPARSE_SOURCE}/LICENSE.txt" "${PREFIX}/share/licenses/suitesparse/"
cp "${SUITESPARSE_SOURCE}/CONTRIBUTOR-LICENSE.txt" "${PREFIX}/share/licenses/suitesparse/"
for package in AMD BTF COLAMD KLU; do
	cp "${SUITESPARSE_SOURCE}/${package}/Doc/License.txt" "${PREFIX}/share/licenses/suitesparse/${package}-License.txt"
done

log "Patching SUNDIALS ${SUNDIALS_VERSION} for CMake 4 compatibility"
"${PYTHON}" - "${SUNDIALS_SOURCE}" <<'PY'
from pathlib import Path
import re
import sys

source = Path(sys.argv[1])
pattern = re.compile(r"(CMAKE_MINIMUM_REQUIRED|cmake_minimum_required)\(VERSION (2\.4|3\.0\.2|3\.1\.3)\)")
patched = []
for path in source.rglob("*.cmake"):
    text = path.read_text()
    updated = pattern.sub(lambda match: f"{match.group(1)}(VERSION 3.5)", text)
    if updated != text:
        path.write_text(updated)
        patched.append(path.relative_to(source))

if not patched:
    raise SystemExit("No SUNDIALS CMake compatibility declarations were patched")

for path in patched:
    print(f"patched {path}")
PY

# Ensure the policy floor also reaches the KLU try_compile subproject, whose
# generated CMakeLists declares an ancient minimum version.
export CMAKE_POLICY_VERSION_MINIMUM=3.5

log "Building SUNDIALS static libraries with MSVC"
SUNDIALS_BUILD="${WORKDIR}/sundials-${SUNDIALS_VERSION}-build"
rm -rf "${SUNDIALS_BUILD}"
cmake -S "${SUNDIALS_SOURCE}" -B "${SUNDIALS_BUILD}" \
	"${CMAKE_GENERATOR_ARGS[@]}" \
	-DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
	-DCMAKE_INSTALL_PREFIX="${PREFIX}" \
	-DCMAKE_INSTALL_LIBDIR=lib \
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
	-DKLU_LIBRARY="${PREFIX}/lib/klu.lib" \
	-DSUNDIALS_INDEX_SIZE=64
cmake --build "${SUNDIALS_BUILD}" --config Release -j "${JOBS}"
cmake --install "${SUNDIALS_BUILD}" --config Release

cp "${SUNDIALS_SOURCE}/LICENSE" "${PREFIX}/share/licenses/sundials/"
cp "${SUNDIALS_SOURCE}/NOTICE" "${PREFIX}/share/licenses/sundials/"

cat >"${PREFIX}/share/cplantbox-native-deps-provenance.txt" <<EOF
CPlantBox source-built native dependencies (Windows/MSVC)

SUNDIALS
  version: ${SUNDIALS_VERSION}
  source: ${SUNDIALS_URL}
  sha256: ${SUNDIALS_SHA256}
  build: static, KLU enabled, SUNDIALS_INDEX_SIZE=64

SuiteSparse
  version: ${SUITESPARSE_VERSION}
  source: ${SUITESPARSE_URL}
  commit: ${SUITESPARSE_COMMIT}
  build: static MSVC via scripts/deps/suitesparse-msvc (SuiteSparse_config, AMD, BTF, COLAMD, KLU)
EOF

log "Installed native dependencies into ${PREFIX}"
