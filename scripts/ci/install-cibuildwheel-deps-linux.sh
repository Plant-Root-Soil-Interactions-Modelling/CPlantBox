#!/usr/bin/env bash
set -euo pipefail

# Install the small set of system tools needed by CPlantBox's manylinux
# cibuildwheel build. The manylinux images already provide most compiler
# tooling, so keep this intentionally minimal.

cmake_is_usable() {
  command -v cmake >/dev/null 2>&1 || return 1
  local version major minor
  version=$(cmake --version | awk '/cmake version/ {print $3; exit}')
  major=${version%%.*}
  minor=${version#*.}
  minor=${minor%%.*}
  [[ ${major:-0} -gt 3 || ( ${major:-0} -eq 3 && ${minor:-0} -ge 5 ) ]]
}

install_with_dnf() {
  local packages=(git curl tar gzip make patch which gcc gcc-c++)
  if ! cmake_is_usable; then
    packages+=(cmake)
  fi
  dnf install -y --setopt=install_weak_deps=False "${packages[@]}"
}

install_with_yum() {
  local packages=(git curl tar gzip make patch which gcc gcc-c++)
  if ! cmake_is_usable; then
    packages+=(cmake)
  fi
  yum install -y "${packages[@]}"
}

install_with_apt() {
  local packages=(git curl tar gzip make patch gcc g++)
  if ! cmake_is_usable; then
    packages+=(cmake)
  fi
  apt-get update
  apt-get install -y --no-install-recommends "${packages[@]}"
}

if command -v dnf >/dev/null 2>&1; then
  install_with_dnf
elif command -v yum >/dev/null 2>&1; then
  install_with_yum
elif command -v apt-get >/dev/null 2>&1; then
  install_with_apt
else
  echo "No supported package manager found for cibuildwheel dependency setup" >&2
  exit 1
fi

cmake --version
