#!/usr/bin/env bash
set -euo pipefail

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "source-native-deps-macos.sh must run on macOS/Darwin" >&2
  exit 2
fi

export MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-11.0}
case " ${CFLAGS:-} " in
  *" -mmacosx-version-min="*) ;;
  *) export CFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} ${CFLAGS:-}" ;;
esac
case " ${CXXFLAGS:-} " in
  *" -mmacosx-version-min="*) ;;
  *) export CXXFLAGS="-mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET} ${CXXFLAGS:-}" ;;
esac

# Reuse the shared SuiteSparse/SUNDIALS source build recipe after applying the
# macOS deployment target flags required for portable macOS wheels.
exec "$(dirname "$0")/source-native-deps.sh" "$@"
