#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

source "$REPO_ROOT/scripts/common/env.sh"
source "$REPO_ROOT/scripts/common/logging.sh"

log_info "[ubuntu] Building test image..."
"$REPO_ROOT/scripts/ubuntu/build_image.sh"

log_info "[ubuntu] Running wheel build + smoke..."
"$REPO_ROOT/scripts/ubuntu/smoke_wheel.sh"

log_info "[ubuntu] Completed."


