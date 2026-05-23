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

Run the standalone smoke helper on a macOS host:

```bash
scripts/wheel/smoke-test-macos.sh
```

To select a specific interpreter, set `CPB_PYTHON`:

```bash
CPB_PYTHON="$(uv python find 3.12)" scripts/wheel/smoke-test-macos.sh
```

The helper:

- builds source native dependencies
- builds a wheel
- audits dependencies with `delocate-listdeps` and `otool`
- installs the wheel into a clean virtual environment
- runs `scripts/wheel/smoke_test.py` outside the source tree

Local `cibuildwheel --platform macos` may require official python.org framework Python package receipts. GitHub-hosted macOS runners are the intended cibuildwheel validation path.
