# macOS wheel builds

CPlantBox macOS wheels are built with the same scikit-build-core and cibuildwheel path used by Linux wheels.

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

The helper is useful for local host validation when `act` cannot execute a real macOS runner. Local `cibuildwheel --platform macos` requires official python.org framework Python receipts; GitHub-hosted macOS runners can install/use those in CI.

The helper:

- builds source native dependencies
- builds a wheel
- audits dependencies with `delocate-listdeps` and `otool`
- installs the wheel into a clean virtual environment
- runs `scripts/wheel/smoke_test.py` outside the source tree

## GitHub Actions workflow

The workflow is in:

```text
.github/workflows/macos-wheels.yml
```

It is written with GitHub Actions-compatible action syntax so the mirrored GitHub repository can run it on GitHub-hosted macOS runners. The job installs a pinned `uv` with `astral-sh/setup-uv`, then uses a `uv`-managed Python environment to run `cibuildwheel` and the wheel audit tools. Local `ACT=true` validation uses `uv` Python interpreters to build, audit, install, and smoke-test all four configured macOS Python wheels on the host.

Forgejo does not directly provide GitHub-hosted macOS runners; if Forgejo needs to display GitHub macOS CI results later, add a status/comment bridge rather than trying to share runners.
