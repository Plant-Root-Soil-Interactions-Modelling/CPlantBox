#!/usr/bin/env bash
set -euo pipefail

# Aggregate Ubuntu platform checks: source-tree smoke and wheel smoke

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

echo "[ubuntu] Running source-tree smoke (dockerized)..."
./scripts/run_ubuntu_tests.sh

echo "[ubuntu] Running wheel build + smoke (dockerized)..."
./scripts/test_wheel_ubuntu.sh

echo "[ubuntu] All Ubuntu checks completed."


