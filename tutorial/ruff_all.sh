#!/usr/bin/env bash

set -euo pipefail

# Find all directories starting with "chapter" at any depth
find . -type d -name "chapter*" | while read -r folder; do
  echo "Processing $folder"

  ruff check --select I "$folder" --fix
  ruff format "$folder"
done

