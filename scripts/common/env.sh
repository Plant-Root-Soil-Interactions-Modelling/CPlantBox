#!/usr/bin/env bash

# Deterministic defaults for local tests/builds
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

# Resolve repo root (call as: REPO_ROOT=$(repo_root))
repo_root() {
  local d
  d="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
  echo "$d"
}


