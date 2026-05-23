# macOS wheel builds

CPlantBox macOS wheels are built by `.github/workflows/wheels.yml` with the same `scikit-build-core` and `cibuildwheel` path used for Linux wheels.

The macOS job builds the runner's native architecture for:

```text
cp311 cp312 cp313 cp314
```

It does not build universal2 wheels or cross-build another macOS architecture.

## Native dependencies

macOS wheels use pinned SuiteSparse/SUNDIALS source builds via:

```bash
scripts/deps/source-native-deps.sh /tmp/cplantbox-native-deps
```

The script builds static archives into a dedicated `cplantbox-*` prefix and avoids Homebrew shared-library runtime dependencies. The wheel config sets the macOS deployment-target flags before invoking it.

## Validation

GitHub-hosted macOS runners are the intended validation path. Local `cibuildwheel --platform macos` may require Python framework installations that are not present on a normal developer machine.

For manual installed-wheel checks, download workflow artifacts and install from a local wheelhouse in a separate test project:

```bash
scripts/download-latest-wheel-artifacts.sh wheelhouse
mkdir -p /tmp/cplantbox-wheel-check
cd /tmp/cplantbox-wheel-check
uv init --bare
uv add --no-index --find-links /path/to/wheelhouse cplantbox
uv run python /path/to/CPlantBox/scripts/wheel/smoke_test.py
```
