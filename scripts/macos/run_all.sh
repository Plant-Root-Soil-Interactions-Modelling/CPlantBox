#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

source "$REPO_ROOT/scripts/common/env.sh"
source "$REPO_ROOT/scripts/common/logging.sh"

log_info "[macOS] Building wheel..."
"$REPO_ROOT/scripts/macos/build/build-wheel.sh"

log_info "[macOS] Testing built wheel..."
"$REPO_ROOT/scripts/macos/test/test-wheel.sh"

log_info "[macOS] Completed build and test."


