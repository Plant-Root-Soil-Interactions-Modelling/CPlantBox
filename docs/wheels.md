# Python wheel builds

CPlantBox platform wheels are built with `scikit-build-core` and `cibuildwheel`.

## Wheel matrix

| Platform | Python tags | Artifact names |
| --- | --- | --- |
| Linux x86_64 manylinux | `cp311`, `cp312`, `cp313`, `cp314` | `cplantbox-manylinux-x86_64-python-3.11` through `cplantbox-manylinux-x86_64-python-3.14` |
| macOS native `macos-latest` architecture | `cp311`, `cp312`, `cp313`, `cp314` | `cplantbox-macos-native-<runner.arch>-python-3.11` through `cplantbox-macos-native-<runner.arch>-python-3.14` |

The macOS job builds the native architecture provided by the GitHub-hosted `macos-latest` runner. It does not cross-build other macOS architectures or universal2 wheels.

Project metadata declares Python 3.11+ plus the core helper-module dependencies (`numpy`, `scipy`, and `matplotlib`). Heavier helper stacks are optional extras: `cplantbox[vis]`, `cplantbox[fipy]`, `cplantbox[mpi]`, or `cplantbox[all]`. ParaView-specific helpers require a ParaView Python environment.

## GitHub workflow

The wheel workflow is:

```text
.github/workflows/wheels.yml
```

It runs on pull requests, pushes to `master`, and manual dispatch. The workflow:

1. builds pinned SuiteSparse/SUNDIALS dependencies from source
2. builds wheels with `cibuildwheel==3.4.1`
3. runs `scripts/wheel/smoke_test.py` against each installed wheel
4. audits wheel portability
5. uploads one artifact per wheel

## Downloading artifacts

Download wheels from the latest completed successful `Wheels` run:

```bash
scripts/download-latest-wheel-artifacts.sh wheelhouse
```

Then install from the local wheelhouse in a separate `uv` project:

```bash
uv init --no-workspace
uv add --find-links /path/to/wheelhouse cplantbox
```

Use `--no-workspace` when creating a temporary test project inside this repository; otherwise uv may add it to the repository workspace and resolve the local source tree instead of the wheel. Do not use `--no-index` unless the wheelhouse also contains all transitive dependencies such as `numpy`, `scipy`, and `matplotlib`.

A later `uv sync` works when both `uv.lock` and the local wheelhouse are present. A `pyproject.toml` alone is not enough until `cplantbox` is published on a package index.

The downloader defaults to the upstream `master` workflow. Override repo/branch when testing a fork:

```bash
CPB_GITHUB_REPO=georgiansarghi/CPlantBox \
CPB_GITHUB_BRANCH=wheels/12-global-polish \
scripts/download-latest-wheel-artifacts.sh wheelhouse
```

## Local cibuildwheel build

Linux x86_64 local builds require Docker:

```bash
python3 -m venv /tmp/cplantbox-cibuildwheel-venv
source /tmp/cplantbox-cibuildwheel-venv/bin/activate
python -m pip install cibuildwheel==3.4.1 auditwheel
python -m cibuildwheel --platform linux --output-dir wheelhouse
CPB_EXPECTED_WHEEL_COUNT=4 AUDITWHEEL=/tmp/cplantbox-cibuildwheel-venv/bin/auditwheel \
  scripts/ci/audit-wheelhouse-linux.sh wheelhouse
```

macOS cibuildwheel validation is intended for GitHub-hosted macOS runners. Local macOS hosts may be unable to run `cibuildwheel --platform macos` unless they have the Python framework installations expected by cibuildwheel.

## Wheel checks

The installed-wheel smoke test verifies that:

- `import plantbox` resolves to the installed wheel, not the source tree
- the compiled extension is imported from site-packages
- `plantbox.data_path()` points to installed model data
- `modelparameter/structural/plant/fspm2023.xml` is present
- a minimal non-graphical plant simulation runs

The audit scripts fail on common portability blockers:

- runtime `libpython` / `Python.framework` dependencies
- absolute runtime search paths
- Homebrew, `/tmp`, user-home, or CI build-directory runtime paths
- unexpected dynamic SuiteSparse/SUNDIALS dependencies

## Publishing status

This PR builds and uploads CI wheel artifacts only. It does not publish to TestPyPI or PyPI.

Wheel artifacts include the source-built SuiteSparse/SUNDIALS license and provenance files under `.dist-info/licenses/third_party/native-dependencies`. Public publishing to TestPyPI/PyPI remains out of scope for this PR.
