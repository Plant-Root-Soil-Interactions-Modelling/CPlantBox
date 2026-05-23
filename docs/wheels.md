# Python wheel builds

CPlantBox platform wheels are built with `scikit-build-core` and `cibuildwheel`.

## Supported wheel matrix

| Platform | Python tags | Artifact names |
| --- | --- | --- |
| Linux x86_64 manylinux | `cp311`, `cp312`, `cp313`, `cp314` | `cplantbox-manylinux-x86_64-python-3.11` through `cplantbox-manylinux-x86_64-python-3.14` |
| macOS native runner architecture | `cp311`, `cp312`, `cp313`, `cp314` | `cplantbox-macos-native-<arch>-python-3.11` through `cplantbox-macos-native-<arch>-python-3.14` |

Linux aarch64, macOS universal2, and Windows wheels are not enabled in this PR.

Project metadata declares:

```toml
requires-python = ">=3.11"
```

## GitHub workflow

The wheel workflow is:

```text
.github/workflows/wheels.yml
```

It builds wheel artifacts on pull requests and manual dispatch:

```bash
gh workflow run wheels.yml --ref <branch>
```

The workflow:

1. builds pinned source native dependencies for each platform
2. builds wheels with `cibuildwheel==3.4.1`
3. runs `scripts/wheel/smoke_test.py` against each installed wheel
4. audits native dependencies
5. uploads one clearly named artifact per wheel

## Local smoke builds

Linux:

```bash
scripts/wheel/smoke-test-linux.sh
```

macOS:

```bash
scripts/wheel/smoke-test-macos.sh
```

The smoke helpers build source native dependencies, build one wheel with the active Python interpreter, install that wheel into a clean virtual environment, and run the installed-package smoke test from outside the source tree.

On macOS, a specific interpreter can be selected with:

```bash
CPB_PYTHON="$(uv python find 3.12)" scripts/wheel/smoke-test-macos.sh
```

## Local cibuildwheel build

Linux x86_64 requires Docker:

```bash
python3 -m venv /tmp/cplantbox-cibuildwheel-venv
source /tmp/cplantbox-cibuildwheel-venv/bin/activate
python -m pip install cibuildwheel==3.4.1 auditwheel
python -m cibuildwheel --platform linux --output-dir wheelhouse
CPB_EXPECTED_WHEEL_COUNT=4 AUDITWHEEL=/tmp/cplantbox-cibuildwheel-venv/bin/auditwheel \
  scripts/ci/audit-wheelhouse-linux.sh wheelhouse
```

macOS cibuildwheel builds are intended for GitHub-hosted macOS runners. Local macOS hosts that do not have official python.org framework Python package receipts should use `scripts/wheel/smoke-test-macos.sh` for host validation.

## Wheel audits

Linux:

```bash
scripts/ci/audit-wheelhouse-linux.sh wheelhouse
```

macOS:

```bash
scripts/ci/audit-wheelhouse-macos.sh wheelhouse
```

The audits fail on common release blockers such as:

- runtime `libpython` / `Python.framework` dependencies
- absolute runtime search paths
- Homebrew, `/tmp`, user-home, or CI build-directory runtime paths
- unexpected dynamic SuiteSparse/SUNDIALS dependencies

## Installed-wheel smoke test

The shared smoke test is:

```text
scripts/wheel/smoke_test.py
```

It verifies that:

- `import plantbox` resolves to the installed wheel, not the source tree
- the compiled extension is imported from site-packages
- `plantbox.data_path()` points to installed package data
- `modelparameter/structural/plant/fspm2023.xml` is present
- a minimal non-graphical plant simulation runs

## Publishing status

This PR builds and uploads CI wheel artifacts only. It does not publish to TestPyPI or PyPI.

Public publishing should wait until third-party license/provenance metadata for bundled/source-built SuiteSparse and SUNDIALS artifacts is complete.
