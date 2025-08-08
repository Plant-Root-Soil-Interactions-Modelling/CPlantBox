#!/usr/bin/env bash
set -euo pipefail

# Build a local wheel using scikit-build-core

python3 -m pip install -U pip build
python3 -m build --wheel
ls -lh dist/*.whl

