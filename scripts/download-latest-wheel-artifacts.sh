#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: scripts/download-latest-wheel-artifacts.sh [OUT_DIR]

Download wheels from the latest completed successful GitHub Actions wheel run and
flatten all .whl files into OUT_DIR.

This intentionally ignores queued/running/failed workflow runs by querying only
runs with status=success.

Requirements:
  gh                  GitHub CLI, installed and authenticated with access to repo artifacts

Environment overrides:
  CPB_GITHUB_REPO     GitHub repo [default: Plant-Root-Soil-Interactions-Modelling/CPlantBox]
  CPB_GITHUB_WORKFLOW Workflow name or file [default: Wheels]
  CPB_GITHUB_BRANCH   Branch to search [default: master]
  CPB_ARTIFACT_NAME   Optional artifact name to download instead of all artifacts

OUT_DIR is deleted and recreated. For safety its basename must be `wheelhouse`
or start with `cplantbox-`.

Examples:
  scripts/download-latest-wheel-artifacts.sh wheelhouse
  CPB_ARTIFACT_NAME=cplantbox-macos-native-ARM64-python-3.11 \
    scripts/download-latest-wheel-artifacts.sh wheelhouse

Install from the downloaded wheelhouse with uv:
  uv add --no-index --find-links wheelhouse cplantbox
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

out_dir=${1:-wheelhouse}
repo=${CPB_GITHUB_REPO:-Plant-Root-Soil-Interactions-Modelling/CPlantBox}
workflow=${CPB_GITHUB_WORKFLOW:-Wheels}
branch=${CPB_GITHUB_BRANCH:-master}
artifact_name=${CPB_ARTIFACT_NAME:-}

command -v python3 >/dev/null 2>&1 || {
  echo "error: python3 is required" >&2
  exit 2
}

abs_path() {
  python3 -c 'import pathlib, sys; print(pathlib.Path(sys.argv[1]).expanduser().resolve(strict=False))' "$1"
}

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(abs_path "${script_dir}/..")
out_dir_abs=$(abs_path "${out_dir}")
out_dir_base=$(basename "${out_dir_abs}")
home_abs=${HOME:+$(abs_path "${HOME}")}

case "${out_dir}" in
  ..|../*|*/..|*/../*)
    echo "error: refusing OUT_DIR containing '..': ${out_dir}" >&2
    exit 2
    ;;
esac

case "${out_dir_base}" in
  wheelhouse|cplantbox-*) ;;
  *)
    echo "error: OUT_DIR basename must be 'wheelhouse' or start with 'cplantbox-': ${out_dir}" >&2
    exit 2
    ;;
esac

case "${out_dir_abs}" in
  /|"${repo_root}"|"${home_abs:-__no_home__}")
    echo "error: refusing to delete unsafe OUT_DIR: ${out_dir_abs}" >&2
    exit 2
    ;;
  "${repo_root}/.."|"${repo_root}/../"*)
    echo "error: refusing to delete path outside the repository via '..': ${out_dir}" >&2
    exit 2
    ;;
esac

command -v gh >/dev/null 2>&1 || {
  cat >&2 <<'EOF'
error: GitHub CLI 'gh' is required.
Install it from https://cli.github.com/ and authenticate with:
  gh auth login
EOF
  exit 2
}

gh auth status --hostname github.com >/dev/null 2>&1 || {
  cat >&2 <<'EOF'
error: GitHub CLI is not authenticated for github.com.
Authenticate with:
  gh auth login
EOF
  exit 2
}

run_id=$(
  gh run list \
    --repo "${repo}" \
    --workflow "${workflow}" \
    --branch "${branch}" \
    --status success \
    --limit 1 \
    --json databaseId \
    --jq '.[0].databaseId'
)

if [[ -z "${run_id}" || "${run_id}" == "null" ]]; then
  echo "error: no successful workflow run found for ${repo}/${workflow} on ${branch}" >&2
  exit 1
fi

rm -rf "${out_dir_abs}"
mkdir -p "${out_dir_abs}"

if [[ -n "${artifact_name}" ]]; then
  echo "Downloading artifact ${artifact_name} from latest successful run ${run_id}"
  gh run download "${run_id}" \
    --repo "${repo}" \
    --name "${artifact_name}" \
    --dir "${out_dir_abs}"
else
  echo "Downloading all artifacts from latest successful run ${run_id}"
  gh run download "${run_id}" \
    --repo "${repo}" \
    --dir "${out_dir_abs}"
fi

# gh stores each artifact in a subdirectory when downloading multiple artifacts.
# Flatten wheels into OUT_DIR for --find-links consumers.
while IFS= read -r -d '' wheel; do
  if [[ "$(dirname "${wheel}")" != "${out_dir_abs}" ]]; then
    mv "${wheel}" "${out_dir_abs}/"
  fi
done < <(find "${out_dir_abs}" -mindepth 2 -type f -name '*.whl' -print0)
find "${out_dir_abs}" -mindepth 1 -type d -empty -delete

wheel_count=$(find "${out_dir_abs}" -maxdepth 1 -type f -name '*.whl' | wc -l | tr -d ' ')
if [[ "${wheel_count}" == "0" ]]; then
  echo "error: no wheels downloaded into ${out_dir_abs}" >&2
  exit 1
fi

find "${out_dir_abs}" -maxdepth 1 -type f -name '*.whl' -print | sort
