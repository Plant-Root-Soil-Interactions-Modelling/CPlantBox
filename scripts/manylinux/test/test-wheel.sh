#!/usr/bin/env bash
set -euo pipefail

# Test a manylinux wheel by running the golden smoke in a matching manylinux test image.
# Usage: scripts/manylinux/test/test-wheel.sh [x86_64|aarch64] [WHEEL_PATH]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$REPO_ROOT"

ARCH="${1:-x86_64}"
WHEEL="${2:-}"
PYTAG_FILTER="${PYTAG_FILTER:-cp310}"

if [[ -z "$WHEEL" ]]; then
  DIR="wheelhouse/linux/manylinux2014_${ARCH}"
  WHEEL=$(ls -1t "$DIR"/*"${PYTAG_FILTER}"*.whl | head -n1)
fi
[[ -f "$WHEEL" ]] || { echo "Wheel not found: $WHEEL" >&2; exit 2; }

# Reuse Ubuntu headless runner for reliable golden rendering
case "$ARCH" in
  x86_64) PLATFORM="linux/amd64";;
  aarch64) PLATFORM="linux/arm64";;
  *) echo "Unsupported arch: $ARCH" >&2; exit 2;;
esac

exec "$REPO_ROOT/scripts/ubuntu/test/test-wheel.sh" "$PLATFORM" "$WHEEL"


