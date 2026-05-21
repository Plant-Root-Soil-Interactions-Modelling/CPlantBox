#!/usr/bin/env bash
set -euo pipefail

wheelhouse=${1:-wheelhouse}
expected_count=${CPB_EXPECTED_WHEEL_COUNT:-}
delocate_listdeps=${DELOCATE_LISTDEPS:-delocate-listdeps}

die() {
  echo "error: $*" >&2
  exit 1
}

[[ "$(uname -s)" == "Darwin" ]] || die "macOS wheel audit must run on Darwin"
command -v "${delocate_listdeps}" >/dev/null 2>&1 || die "delocate-listdeps not found: ${delocate_listdeps}"
command -v python3 >/dev/null 2>&1 || die "python3 not found"

[[ -d "${wheelhouse}" ]] || die "wheelhouse does not exist: ${wheelhouse}"

wheels=()
while IFS= read -r wheel; do
  wheels+=("${wheel}")
done < <(find "${wheelhouse}" -maxdepth 1 -type f -name '*.whl' | sort)
wheel_count=${#wheels[@]}
[[ "${wheel_count}" -gt 0 ]] || die "no wheels found in ${wheelhouse}"

if [[ -n "${expected_count}" && "${wheel_count}" != "${expected_count}" ]]; then
  die "expected ${expected_count} wheels, found ${wheel_count}"
fi

tmp_root=$(mktemp -d "${TMPDIR:-/tmp}/cplantbox-macos-wheel-audit.XXXXXXXX")
cleanup() {
  case "${tmp_root}" in
    /tmp/cplantbox-macos-wheel-audit.*|${TMPDIR:-/tmp}/cplantbox-macos-wheel-audit.*) rm -rf "${tmp_root}" ;;
    *) echo "warning: refusing to remove unexpected temp path: ${tmp_root}" >&2 ;;
  esac
}
trap cleanup EXIT

for wheel in "${wheels[@]}"; do
  wheel_name=$(basename "${wheel}")
  echo "::group::audit ${wheel_name}"

  case "${wheel_name}" in
    *macosx*.whl) ;;
    *) die "wheel does not carry an expected macOS tag: ${wheel_name}" ;;
  esac

  deps_log="${tmp_root}/${wheel_name}.delocate-listdeps.txt"
  "${delocate_listdeps}" "${wheel}" | tee "${deps_log}"
  ! grep -E '(^|[[:space:]])(/opt/homebrew|/usr/local/Cellar|/tmp|/private/tmp|/Users)/' "${deps_log}" || die "delocate found a non-portable dependency path in ${wheel_name}"

  extract_dir="${tmp_root}/${wheel_name}.contents"
  mkdir -p "${extract_dir}"
  python3 - "${wheel}" "${extract_dir}" <<'PY'
from pathlib import Path
import sys
import zipfile

wheel = Path(sys.argv[1])
dest = Path(sys.argv[2])
with zipfile.ZipFile(wheel) as archive:
    archive.extractall(dest)
PY

  mach_objects=()
  while IFS= read -r object; do
    mach_objects+=("${object}")
  done < <(find "${extract_dir}" -type f \( -name '*.so' -o -name '*.dylib' \) | sort)
  [[ "${#mach_objects[@]}" -gt 0 ]] || die "no Mach-O extension/library files found in ${wheel_name}"

  otool_log="${tmp_root}/${wheel_name}.otool.txt"
  rpath_log="${tmp_root}/${wheel_name}.rpaths.txt"
  : > "${otool_log}"
  : > "${rpath_log}"
  for object in "${mach_objects[@]}"; do
    echo "### ${object#${extract_dir}/}" | tee -a "${otool_log}"
    otool -L "${object}" | tee -a "${otool_log}"
    echo "### ${object#${extract_dir}/}" | tee -a "${rpath_log}"
    otool -l "${object}" | awk '/cmd LC_RPATH/{show=1; next} show && /path /{print; show=0}' | tee -a "${rpath_log}"
  done

  ! grep -E '(^|[[:space:]])(/opt/homebrew|/usr/local/Cellar|/tmp|/private/tmp|/Users)/' "${otool_log}" || die "otool found a non-portable dependency path in ${wheel_name}"
  ! grep -E '(^|[[:space:]])path (/opt/homebrew|/usr/local/Cellar|/tmp|/private/tmp|/Users)/' "${rpath_log}" || die "otool found a non-portable runtime search path in ${wheel_name}"
  ! grep -E 'Python\.framework|libpython[0-9]+\.[0-9]+' "${otool_log}" || die "otool found a Python library dependency in ${wheel_name}"

  echo "::endgroup::"
done

echo "Audited ${wheel_count} macOS wheel(s) in ${wheelhouse}"
