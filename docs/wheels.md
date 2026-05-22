# Python wheel support

This document describes the wheel build and release-readiness plan for CPlantBox. It complements the source-build instructions in `README.md` and the dependency policy in `docs/native-dependencies.md`.

## Current status

Wheel support is being validated in CI before public publishing. Do not treat wheel artifacts from CI as an official release unless they are attached to an approved release and the release checklist in `docs/release.md` is complete.

Publishing release wheels remains blocked until third-party license and provenance metadata is complete for bundled/source-built SuiteSparse and SUNDIALS artifacts.

## Supported Python versions

The project metadata declares:

```toml
requires-python = ">=3.11"
```

The intended wheel matrix is:

| Platform | Python tags | Status |
| --- | --- | --- |
| Linux x86_64 manylinux | `cp311`, `cp312`, `cp313`, `cp314` | Built and audited in Forgejo/Gitea CI |
| Linux aarch64 manylinux | `cp311`, `cp312`, `cp313`, `cp314` | Deferred until native ARM64 runner capacity or acceptable QEMU runtime is available |
| macOS native runner architecture | `cp311` | Build workflow added for future GitHub-hosted macOS runners; locally validated with `act` host execution |
| macOS universal2 | not enabled | Deferred |
| Windows | not available | Out of scope |

## Installing wheels

After official wheels are published, users should be able to install CPlantBox with:

```bash
python -m pip install cplantbox
```

For TestPyPI validation, maintainers can use:

```bash
python -m pip install \
  --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple/ \
  cplantbox
```

If no compatible wheel exists for the user's platform or Python version, `pip` may try a source build or report that no matching distribution is available. Use Python 3.11 or newer and one of the supported wheel platforms above.

## Local wheel build smoke tests

The preferred local validation path is the platform helper script. These helpers build pinned source native dependencies, build a wheel, install it into a clean virtual environment, and run the installed-package smoke test from outside the source tree.

Linux:

```bash
scripts/wheel/smoke-test-linux.sh
```

macOS:

```bash
scripts/wheel/smoke-test-macos.sh
```

Both helpers accept environment overrides such as `CPB_WHEELHOUSE`, `CPB_NATIVE_DEPS_PREFIX`, and `CPB_DEPS_JOBS`. Use dedicated `cplantbox-*` paths because the helpers delete and recreate their output directories.

## Local cibuildwheel commands

Linux x86_64 wheels can be built with Docker and cibuildwheel:

```bash
python3 -m venv /tmp/cplantbox-cibuildwheel-venv
/tmp/cplantbox-cibuildwheel-venv/bin/python -m pip install cibuildwheel==3.4.1 auditwheel
/tmp/cplantbox-cibuildwheel-venv/bin/python -m cibuildwheel --platform linux --output-dir wheelhouse
CPB_EXPECTED_WHEEL_COUNT=4 \
  AUDITWHEEL=/tmp/cplantbox-cibuildwheel-venv/bin/auditwheel \
  scripts/ci/audit-wheelhouse-linux.sh wheelhouse
```

On macOS, local `cibuildwheel --platform macos` expects official python.org framework Python package receipts. If the host does not have those receipts, use the macOS helper or the workflow's `ACT=true` host-validation path instead:

```bash
act workflow_dispatch \
  -W .github/workflows/macos-wheels.yml \
  -j macos-wheel \
  -P macos-latest=-self-hosted \
  --env ACT=true
```

## CI wheel artifacts

Current artifact names:

| Workflow | Artifact | Contents |
| --- | --- | --- |
| `.gitea/workflows/ubuntu-baseline-build.yml` | `cplantbox-manylinux-x86_64` | Linux x86_64 manylinux wheels for the configured Python matrix |
| `.github/workflows/macos-wheels.yml` | `cplantbox-macos-native` | macOS native-architecture wheel |
| `.github/workflows/release-wheels.yml` | `cplantbox-release-readiness` | Dry-run release readiness summary, not wheels |

Generated wheels belong in CI artifacts or local `wheelhouse/` directories. They must not be committed to git.

## Inspecting wheels

Linux:

```bash
auditwheel show wheelhouse/*.whl
scripts/ci/audit-wheelhouse-linux.sh wheelhouse
```

macOS:

```bash
delocate-listdeps wheelhouse/*.whl
scripts/ci/audit-wheelhouse-macos.sh wheelhouse
```

The audits check for common release blockers, including `libpython` dependencies, host-specific runtime paths, and accidental dynamic SuiteSparse/SUNDIALS dependencies.

## Installed-wheel smoke test expectations

The shared smoke test is `scripts/wheel/smoke_test.py`. It verifies that:

- `import plantbox` resolves to the installed wheel, not the source tree
- the compiled extension is imported from site-packages
- `plantbox.data_path()` points to installed package data
- `modelparameter/structural/plant/fspm2023.xml` is present
- a minimal non-graphical plant simulation runs

## Troubleshooting

### Import uses the source tree

Run the smoke test from outside the repository, for example from `/tmp`. The smoke test removes CPlantBox source-tree entries from `sys.path`, but running ad-hoc checks from the repository can still hide packaging mistakes.

### Missing native library at import time

Inspect the wheel with `auditwheel show` on Linux or `delocate-listdeps`/`otool -L` on macOS. Release wheels should not depend on arbitrary libraries from Homebrew, `/tmp`, a user home directory, or CI build directories.

### Unexpected `libpython` dependency

The Python extension should link as an extension module, not as an embedded-Python application. Linux and macOS audit helpers fail if a wheel depends on `libpython`/`Python.framework`.

### Missing package data

Check:

```python
import plantbox as pb
print(pb.data_path())
```

The path should exist inside the installed `plantbox` package and contain the `modelparameter` tree.

## Rendering and golden tests

Rendering and golden-image tests are not required wheel gates yet. The current required wheel gate is the installed-package import/data/minimal-simulation smoke test.

Recommended progression:

1. Keep import/package-data/minimal-simulation smoke tests required on every wheel build.
2. Add optional or nightly VTK/headless rendering smoke tests later.
3. Add platform-specific golden image comparisons only after rendering output is stable across platforms.
4. Write generated images to temporary/artifact directories, never into tracked `test/golden` paths during CI.
