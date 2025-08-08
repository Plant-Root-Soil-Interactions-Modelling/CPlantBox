#!/usr/bin/env bash
set -euo pipefail

# Tag a release (annotated), optionally push, and build Ubuntu wheel artifacts.
#
# Usage:
#   scripts/release/release_tag.sh vX.Y[.Z][-rcN] [-m "message"] [--push] [--skip-tests]
#
# Examples:
#   scripts/release/release_tag.sh v2.1 -m "WIP: packaging modernization kickoff"
#   scripts/release/release_tag.sh v2.1.1 --push

if [[ ${1-} == "" ]]; then
  echo "ERROR: Version tag required (e.g., v2.1 or v2.1.0)" >&2
  exit 2
fi

VERSION="$1"; shift || true
MESSAGE="Release ${VERSION}"
PUSH=false
SKIP_TESTS=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    -m|--message)
      MESSAGE="$2"; shift 2 ;;
    --push)
      PUSH=true; shift ;;
    --skip-tests)
      SKIP_TESTS=true; shift ;;
    *)
      echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

if ! [[ "$VERSION" =~ ^v[0-9]+(\.[0-9]+){1,2}(-rc[0-9]+)?$ ]]; then
  echo "ERROR: Version must look like v2.1 or v2.1.0 or v2.1-rc1" >&2
  exit 2
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

echo "Ensuring clean working tree..."
if ! git diff-index --quiet HEAD --; then
  echo "ERROR: Working tree is not clean. Commit or stash changes first." >&2
  exit 3
fi

echo "Fetching tags..."
git fetch --tags --quiet || true

if git rev-parse -q --verify "refs/tags/${VERSION}" >/dev/null; then
  echo "ERROR: Tag ${VERSION} already exists" >&2
  exit 4
fi

if [[ "$SKIP_TESTS" == false ]]; then
  echo "Running aggregated Ubuntu checks before tagging..."
  ./scripts/ubuntu/run_all.sh
fi

echo "Creating annotated tag ${VERSION}..."
git tag -a "$VERSION" -m "$MESSAGE"

if [[ "$PUSH" == true ]]; then
  echo "Pushing tag ${VERSION} to origin..."
  git push origin "$VERSION"
fi

echo "Building Ubuntu wheel artifact for ${VERSION}..."
./scripts/test_wheel_ubuntu.sh

ART_DIR="release_artifacts/${VERSION}/ubuntu"
mkdir -p "$ART_DIR"
cp -v dist/*.whl "$ART_DIR"/

echo "Release ${VERSION} prepared. Artifacts in $ART_DIR"


