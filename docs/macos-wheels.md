# macOS wheel builds

CPlantBox macOS wheels are built by `.github/workflows/wheels.yml` with the same scikit-build-core and cibuildwheel path used for Linux wheels.

The current macOS job builds the Python matrix through Python 3.14 on the runner's native architecture:

- Python tags: `cp311`, `cp312`, `cp313`, `cp314`
- Apple Silicon runner: `arm64`
- Intel runner: `x86_64`

Universal2 wheels and cross-architecture builds are intentionally deferred.

## Native dependencies

macOS wheels use source-built, pinned SuiteSparse/SUNDIALS dependencies via:

```bash
scripts/deps/source-native-deps-macos.sh /tmp/cplantbox-native-deps-macos
```

The script builds static archives into a dedicated `cplantbox-*` prefix. This avoids relying on Homebrew shared libraries at wheel runtime.

## Local validation on macOS

GitHub-hosted macOS runners are the intended cibuildwheel validation path. Local `cibuildwheel --platform macos` may require official python.org framework Python package receipts.

For manual installed-wheel checks, download workflow artifacts and install them from a local wheelhouse:

```bash
scripts/download-latest-wheel-artifacts.sh wheelhouse
uv add --find-links wheelhouse cplantbox
uv run python scripts/wheel/smoke_test.py
```
