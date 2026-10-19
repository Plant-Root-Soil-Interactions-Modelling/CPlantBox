#!/usr/bin/env bash
set -euo pipefail

wheelhouse=${1:-wheelhouse}
expected_count=${CPB_EXPECTED_WHEEL_COUNT:-}
auditwheel=${AUDITWHEEL:-auditwheel}

die() {
  echo "error: $*" >&2
  exit 1
}

command -v "${auditwheel}" >/dev/null 2>&1 || die "auditwheel not found: ${auditwheel}"
command -v readelf >/dev/null 2>&1 || die "readelf not found"
command -v python3 >/dev/null 2>&1 || die "python3 not found"

[[ -d "${wheelhouse}" ]] || die "wheelhouse does not exist: ${wheelhouse}"

mapfile -t wheels < <(find "${wheelhouse}" -maxdepth 1 -type f -name '*.whl' | sort)
wheel_count=${#wheels[@]}
[[ "${wheel_count}" -gt 0 ]] || die "no wheels found in ${wheelhouse}"

if [[ -n "${expected_count}" && "${wheel_count}" != "${expected_count}" ]]; then
  die "expected ${expected_count} wheels, found ${wheel_count}"
fi

tmp_root=$(mktemp -d "${TMPDIR:-/tmp}/cplantbox-wheel-audit.XXXXXXXX")
cleanup() {
  case "${tmp_root}" in
    /tmp/cplantbox-wheel-audit.*|${TMPDIR:-/tmp}/cplantbox-wheel-audit.*) rm -rf "${tmp_root}" ;;
    *) echo "warning: refusing to remove unexpected temp path: ${tmp_root}" >&2 ;;
  esac
}
trap cleanup EXIT

for wheel in "${wheels[@]}"; do
  wheel_name=$(basename "${wheel}")
  echo "::group::audit ${wheel_name}"

  case "${wheel_name}" in
    *manylinux*x86_64.whl) ;;
    *) die "wheel does not carry an expected manylinux x86_64 tag: ${wheel_name}" ;;
  esac

  audit_log="${tmp_root}/${wheel_name}.auditwheel.txt"
  "${auditwheel}" show "${wheel}" | tee "${audit_log}"
  ! grep -E 'libpython[0-9]+\.[0-9]+\.so' "${audit_log}" || die "auditwheel found a libpython dependency in ${wheel_name}"
  ! grep -E '/(home|tmp|var|private|Users)/' "${audit_log}" || die "auditwheel output contains a host-specific path in ${wheel_name}"

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

  python3 - "${extract_dir}" <<'PY'
from pathlib import Path
import sys

root = Path(sys.argv[1])
dist_infos = list(root.glob("*.dist-info"))
if len(dist_infos) != 1:
    raise SystemExit(f"expected one .dist-info directory, found {len(dist_infos)}")
license_root = dist_infos[0] / "licenses" / "third_party" / "native-dependencies"
required = [
    "cplantbox-native-deps-provenance.txt",
    "sundials/LICENSE",
    "sundials/NOTICE",
    "suitesparse/LICENSE.txt",
]
missing = [name for name in required if not (license_root / name).is_file()]
if missing:
    raise SystemExit(f"missing native dependency license/provenance files: {missing}")
PY

  mapfile -t shared_objects < <(find "${extract_dir}" -type f \( -name '*.so' -o -name '*.so.*' \) | sort)
  [[ "${#shared_objects[@]}" -gt 0 ]] || die "no shared objects found in ${wheel_name}"

  readelf_log="${tmp_root}/${wheel_name}.readelf.txt"
  : > "${readelf_log}"
  for object in "${shared_objects[@]}"; do
    echo "### ${object#${extract_dir}/}" | tee -a "${readelf_log}"
    readelf -d "${object}" | tee -a "${readelf_log}"
  done

  ! grep -E 'NEEDED.*libpython[0-9]+\.[0-9]+\.so' "${readelf_log}" || die "readelf found a libpython dependency in ${wheel_name}"
  ! grep -E 'NEEDED.*(libsundials|libklu|libamd|libcolamd|libbtf|libsuitesparseconfig)' "${readelf_log}" || die "readelf found an unexpected dynamic native dependency in ${wheel_name}"
  ! grep -E '(RPATH|RUNPATH).*[[]/' "${readelf_log}" || die "readelf found an absolute runtime search path in ${wheel_name}"
  ! grep -E 'NEEDED.*[[]/' "${readelf_log}" || die "readelf found an absolute NEEDED path in ${wheel_name}"

  echo "::endgroup::"
done

echo "Audited ${wheel_count} wheel(s) in ${wheelhouse}"
