# Wheel release readiness

This document defines the operational release plan for CPlantBox Python wheels.

## Release status

Wheel publishing is not yet enabled. CI can build and audit wheel artifacts, but public TestPyPI/PyPI publishing should remain disabled until all release gates are complete.

Current blocker:

- Complete third-party license and provenance metadata for source-built/bundled SuiteSparse and SUNDIALS artifacts before publishing wheels that contain their compiled code.

## Supported release matrix

Initial release target:

| Platform | Python versions | Release status |
| --- | --- | --- |
| Linux x86_64 manylinux | 3.11, 3.12, 3.13, 3.14 | Ready for CI validation |
| macOS native architecture | 3.11 | Ready for GitHub-hosted macOS validation |
| Linux aarch64 manylinux | 3.11, 3.12, 3.13, 3.14 | Deferred |

The project has no Python upper bound in `pyproject.toml`; CI pins explicit wheel tags so new Python releases are added intentionally.

## Linux aarch64 decision

Linux aarch64 wheels are explicitly deferred for the first release-readiness pass.

Preferred enablement order:

1. Add a native ARM64 Linux runner and build `manylinux_aarch64` wheels without emulation.
2. If no native runner is available, test QEMU-based `cibuildwheel` builds and keep them only if runtime is acceptable.
3. Keep aarch64 deferred rather than blocking x86_64/macOS releases indefinitely.

Before enabling aarch64, validate all of the following:

- source-built SuiteSparse/SUNDIALS completes in the target environment
- `scripts/wheel/smoke_test.py` passes from an installed wheel
- `auditwheel show` and `scripts/ci/audit-wheelhouse-linux.sh` pass
- CI runtime is acceptable for the regular release workflow
- artifacts are named separately from x86_64 wheels, for example `cplantbox-manylinux-aarch64`

## Release workflow stages

Recommended release process:

1. Pull request CI builds and smoke-tests selected wheels.
2. Main branch CI builds and audits the full supported matrix, but does not publish.
3. Tag CI builds and audits the full supported matrix, uploads release artifacts, and optionally publishes to TestPyPI.
4. A protected environment/manual approval publishes to PyPI after artifact and metadata review.

A dry-run release skeleton lives in:

```text
.github/workflows/release-wheels.yml
```

It generates a release-readiness summary artifact and intentionally does not publish packages.

## Release checklist

Before publishing TestPyPI or PyPI wheels:

- [ ] Linux x86_64 wheels build for `cp311`, `cp312`, `cp313`, and `cp314`.
- [ ] macOS wheels build on a GitHub-hosted macOS runner for the configured matrix.
- [ ] Each wheel installs into a clean virtual environment.
- [ ] `scripts/wheel/smoke_test.py` passes from outside the repository tree.
- [ ] `plantbox.data_path()` points to installed package data.
- [ ] Linux wheels pass `auditwheel show` and `scripts/ci/audit-wheelhouse-linux.sh`.
- [ ] macOS wheels pass `delocate-listdeps`/`otool` audit via `scripts/ci/audit-wheelhouse-macos.sh`.
- [ ] No wheel depends on `libpython`/`Python.framework`.
- [ ] No wheel depends on Homebrew, `/tmp`, a user home directory, or CI build-directory runtime paths.
- [ ] Source-built SuiteSparse/SUNDIALS provenance and license metadata is complete.
- [ ] Generated wheels are uploaded as CI/release artifacts and are not committed to git.
- [ ] Release notes state the supported Python/platform matrix and known limitations.

## TestPyPI and PyPI publishing plan

Prefer trusted publishing over long-lived API tokens.

Recommended setup:

1. Create TestPyPI and PyPI projects for `cplantbox`.
2. Configure trusted publishers for the GitHub repository, release workflow, and protected environments.
3. Add a TestPyPI publish job that runs only on tags or manual dispatch after all build/audit jobs pass.
4. Add a PyPI publish job behind a protected environment and manual approval.
5. Keep TestPyPI and PyPI publishing disabled until the license/provenance gate is complete.

Future publishing action shape:

```yaml
- name: Publish to TestPyPI
  uses: pypa/gh-action-pypi-publish@release/v1
  with:
    repository-url: https://test.pypi.org/legacy/
    packages-dir: wheelhouse
```

Do not add this as an active publishing step until trusted publishing and metadata gates are ready.

## Artifact conventions

Use explicit artifact names per platform/architecture:

| Artifact | Contents |
| --- | --- |
| `cplantbox-manylinux-x86_64` | Linux x86_64 wheel matrix |
| `cplantbox-manylinux-aarch64` | Future Linux aarch64 wheel matrix |
| `cplantbox-macos-native` | macOS native-architecture wheels |
| `cplantbox-release-readiness` | Dry-run release summary |

Release artifacts should include only generated distribution files and release summaries. Build directories, virtual environments, `_skbuild/`, and local dependency work directories must not be committed.

## Maintainer commands

Linux local smoke build:

```bash
scripts/wheel/smoke-test-linux.sh
```

macOS local smoke build:

```bash
scripts/wheel/smoke-test-macos.sh
```

Linux cibuildwheel matrix:

```bash
python3 -m venv /tmp/cplantbox-cibuildwheel-venv
/tmp/cplantbox-cibuildwheel-venv/bin/python -m pip install cibuildwheel==3.4.1 auditwheel
/tmp/cplantbox-cibuildwheel-venv/bin/python -m cibuildwheel --platform linux --output-dir wheelhouse
CPB_EXPECTED_WHEEL_COUNT=4 \
  AUDITWHEEL=/tmp/cplantbox-cibuildwheel-venv/bin/auditwheel \
  scripts/ci/audit-wheelhouse-linux.sh wheelhouse
```

macOS workflow validation with local `act` host execution:

```bash
act workflow_dispatch \
  -W .github/workflows/macos-wheels.yml \
  -j macos-wheel \
  -P macos-latest=-self-hosted \
  --env ACT=true
```

## Rendering and golden validation

Rendering/golden tests are not release gates yet.

Policy:

- required for every wheel build: installed-package import, package-data check, minimal simulation smoke test
- optional or nightly: VTK/headless rendering smoke test
- release-gating only after stable: platform-specific golden references
- generated images must go to temporary directories or CI artifacts, never tracked `test/golden` paths

This avoids blocking wheel availability on rendering differences that are not yet understood across Linux/macOS/Python versions.

## Known limitations

- Linux aarch64 is deferred.
- macOS starts with a native `cp311` wheel; universal2 and a larger macOS Python matrix are deferred.
- Windows wheels are not planned in this stack.
- Public publishing is blocked on third-party license/provenance metadata.
