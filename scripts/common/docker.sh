#!/usr/bin/env bash

# Ensure a buildx builder exists and is selected
ensure_buildx_builder() {
  local name=${1:-cplantbox-builder}
  if ! docker buildx inspect "$name" >/dev/null 2>&1; then
    docker buildx create --use --name "$name" >/dev/null
  else
    docker buildx use "$name" >/dev/null
  fi
}

# Build an image with buildx --load to local docker
docker_build_image() {
  local image_name=$1 dockerfile=$2 platform=$3 context=${4:-.}
  docker buildx build --platform "$platform" -f "$dockerfile" -t "$image_name" --load "$context"
}

# Run a command in a container mounting the project root
docker_run_project() {
  local image=$1 platform=$2 workdir=${3:-/src} cmd=$4
  docker run --rm -t --platform "$platform" -v "$(pwd)":/src -w "$workdir" "$image" bash -lc "$cmd"
}


