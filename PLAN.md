## CPlantBox Packaging & Build Modernization Plan

This living document tracks our plan to make CPlantBox easy to install via pip on multiple platforms while keeping researcher workflows intact. We will iterate in small, safe steps and use an Ubuntu amd64 headless smoke test as our guardrail.

### Current state (2025)

- Python API is provided via a pybind11 C++ extension named `plantbox` built by CMake (`src/PyPlantBox.cpp`).
- Heavy native dependencies: SuiteSparse (AMD/COLAMD/BTF/KLU) and SUNDIALS (CVODE and sun* components), plus vendored tinyxml2, gauss_legendre, PiafMunch, Eigen.
- CMake:
  - Old Python detection (`find_package(PythonInterp/PythonLibs 3.7)`).
  - Non-portable flags (e.g., default Debug, and `-march=native` in release flags).
  - Manually forced output directories in places.
- Tutorials/tests import `plantbox` and access resources via relative paths to `modelparameter/**`.
- No modern packaging (`pyproject.toml`) yet. A historic PR attempted `setup.py`-based wheels and Debian packages but is stale and not mergeable as-is.
- CI not configured for multi-platform wheels.

### Our goals

1. Deliver pip-installable wheels for common platforms (Linux manylinux x86_64/aarch64, macOS x86_64/arm64). Windows optional later.
2. Keep researcher UX: fast local builds, tutorials, and examples remain usable; ship `modelparameter/**` in wheels with a stable access helper.
3. Modernize packaging and build with minimal risk: adopt `pyproject.toml` + scikit-build-core and improve CMake hygiene.
4. Establish reliable headless smoke tests to prevent regressions.

### Guardrails (baseline)

- A dedicated Ubuntu amd64 headless smoke test is available:
  - Image: `docker/Dockerfile.ubuntu-test-env` (build toolchain + python deps, headless).
  - Runner: `scripts/run_ubuntu_tests.sh` uses docker buildx (hardcoded `linux/amd64`) to build and run:
    - Curated pytest subset (headless-safe).
    - A small inline example that simulates and writes a VTP file (no VTK rendering).

---

## Phased Plan and Tasks

We track tasks with checkboxes. Check off when done; add notes/links to changes as we progress.

### Phase 0 — Guardrails (done)

- [x] Add Ubuntu 22.04 testing image with toolchain and minimal python deps: `docker/Dockerfile.ubuntu-test-env`.
- Details: Installs build-essential, cmake, git, libgl1, and headless Python deps (pytest, numpy, scipy, matplotlib, pandas, vtk). Keeps image lean and GUI-free.
- Caveats: No GUI stack; VTK used only for file I/O. Heavy GUI tutorials intentionally excluded from smoke.
- [x] Add `scripts/run_ubuntu_tests.sh` to build with buildx (`linux/amd64`) and run headless smoke tests.
- Details: Uses docker buildx and a dedicated builder to build an amd64 image on Apple Silicon, then runs tests with `--platform linux/amd64`.
- Caveats: Requires Docker Desktop with buildx/QEMU. First run may download emulation components and be slower.
- [x] Curate a headless-safe pytest subset and a headless inline simulation (avoid VTK GUI).
- Details: Runs a known-good subset of unit tests and a minimal inline simulation that writes a VTP file, avoiding rendering.
- Caveats: Full test suite still includes GUI/VTK-heavy cases which are not executed in the headless smoke.

Notes: This is our single-command regression check. Keep it green before/after changes.

### Phase 1 — Repo hygiene and small, safe fixes

- [x] Modernize CMake Python detection (use `find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)` or `FindPython3`).
  - Details: Switched to `find_package(Python COMPONENTS Interpreter Development.Module Development.Embed REQUIRED)` and kept vendored pybind11 (`add_subdirectory(external/pybind11)`). `Development.Embed` ensures `Python::Python` exists so `pybind11_add_module` succeeds.
  - Caveats: Requires Python dev headers in the build environment; covered by Docker image via `python3-dev`.
- [x] Default build type to Release for packaging flows; make Debug explicit.
  - Details: Root `CMakeLists.txt` now sets `CMAKE_BUILD_TYPE=Release` if not provided.
  - Caveats: Developers can still pass `-DCMAKE_BUILD_TYPE=Debug` locally.
- [x] Remove CPU-specific flags like `-march=native`; keep portable optimization flags.
  - Details: Removed `-march=native` from release flags; kept `-O3` and safe opts.
  - Caveats: Small performance drop versus native builds; acceptable for portable wheels. Consider adding an opt-in CMake toggle later for local dev.
- [x] Avoid forcing output directories for libs/modules unless required by wheel build tooling.
  - Details: Commented out `LIBRARY_OUTPUT_DIRECTORY` overrides in `src/CMakeLists.txt` for `CPlantBox` and `plantbox`; use CMake defaults for portability.
  - Caveats: Packaging will manage install locations; local builds still place artifacts under build tree as usual.
- [x] Introduce pytest markers to tag VTK-dependent tests; skip automatically in headless environments.
  - Details: Added `pytest.ini` with markers `vtk` and `gui`, and `test/conftest.py` now auto-skips `vtk` tests when VTK is unavailable and `gui` tests when no display is present. Next step (optional) is to annotate specific tests with these markers as needed.
- [x] Add `ruff` (conservative config: E,F,I) and `pre-commit` for python files; document opt-in for devs.
  - Details: Added `ruff.toml` and `.pre-commit-config.yaml` with ruff (E,F,I), black (manual), and basic hooks. Devs can `pip install pre-commit && pre-commit install`.
- [x] Add `clang-format` config (do not enforce globally yet; document usage for contributors).
  - Details: Added `.clang-format` with a lightweight, non-enforced style.
- [x] Disable LTO/IPO to avoid toolchain issues in container builds.
  - Details: Set `CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF` in `src/CMakeLists.txt` to resolve `lto-wrapper`/jobs server errors seen under emulation.
  - Caveats: Minor potential link-time optimization loss; acceptable for reliability in cross-arch container builds.

  New small hygiene tasks (added):
  - [x] Set an explicit minimum CMake version compatible with modern `FindPython` (e.g., `cmake_minimum_required(VERSION 3.18)` or higher) and document it.
    - Details: Updated root `CMakeLists.txt` to `cmake_minimum_required(VERSION 3.18)`.
    - Caveats: Older CMake versions (<3.18) will now fail fast; aligns with modern Python discovery.
  - [x] Add a top-level `.dockerignore` to speed Docker builds and avoid copying large artifacts (e.g., `build/`, `__pycache__/`, `.git/`, `tutorial/results/`).
    - Details: Added `.dockerignore` excluding build outputs, caches, large artifacts, and test outputs.
    - Caveats: If Docker builds need additional files later, update `.dockerignore` accordingly.
  - [x] Make test runs deterministic by default in scripts (export `OMP_NUM_THREADS=1`) and document this for local dev.
    - Details: `scripts/run_ubuntu_tests.sh` and the `CTest` smoke hook set `OMP_NUM_THREADS=1`.
    - Caveats: Users can override locally; ensures stable timings under CI/emulation.
  - [ ] Ensure vendored `external/pybind11` is present and at a known-good tag; add a short note in docs and a friendly check in the build script if missing. (Moved to Phase 2)
  - [x] Add a simple CTest hook that runs the curated pytest subset (`ctest -R smoke`) to make `cmake --build . && ctest` work out-of-the-box.
    - Details: Root `CMakeLists.txt` defines `smoke_pytest` with `python3 -m pytest -q` on a curated set; sets `PYTHONPATH` and `OMP_NUM_THREADS` via test properties.
    - Caveats: Requires Python available in PATH; uses the source tree (not installed wheel) for speed.
  - [x] Ensure tests do not leave artifacts in the repo; add cleanup and ignore rules.
    - Details: Route test-generated files (e.g., `human.xml`, `leaf.xml`, `organ.xml`, `organism.rsml`, `root.xml`, `rs_parameters.xml`, `seed.xml`, `stem.xml`, `test_rootsystem.rsml`) into per-test temp dirs via `tmp_path` and/or chdir fixtures; add autouse teardown to remove legacy outputs; add `.gitignore` patterns as a fallback.
    - Caveats: Some tutorials intentionally write outputs; keep those under `tutorial/results/` (already ignored) and out of repo root.

### Phase 2 — Introduce modern packaging (pyproject + scikit-build-core)

- [x] Add `pyproject.toml` using `scikit-build-core` as the build backend.
  - Details: Added minimal `pyproject.toml` with scikit-build-core; CMake/Ninja injected by backend; Release build via tool config.
  - Caveats: Build is green on Ubuntu container; cross-platform wheels still pending in later phases.
- [x] Configure extension build via CMake without setup.py; verify `pip wheel .` works in the Ubuntu test container.
  - Details: Wheel builds and links successfully under `scikit-build-core`. Implemented out-of-source fallbacks for vendored static libs in `src/CMakeLists.txt`, added `-DSKBUILD=ON` to switch `CPlantBox` to STATIC for wheel builds, and fixed SuiteSparse static link order. Produced wheel (~1.7 MB) installs and runs the headless smoke.
- [x] Keep runtime python deps minimal (e.g., `numpy`), avoid optional heavy ones (`vtk`) in core install.
  - Details: No runtime deps declared yet; tutorials/tests keep optional packages.
- [x] Prove smoke from fresh venv in-container using the built wheel.
  - Details: `scripts/test_wheel_ubuntu.sh` builds the wheel, installs into a temp venv, imports `plantbox`, runs a tiny headless simulation; all green in Ubuntu container.

  Additional Phase 2 tasks (added):
  - [x] Ensure vendored `external/pybind11` is present and at a known-good tag; document the expected version and a friendly check in the wheel build helper.
    - Details: Added a fail-fast check in `src/CMakeLists.txt` that errors with guidance if `src/external/pybind11` is missing; submodule is initialized in scripts.
  - [x] Choose and implement package layout for the Python distribution:
    - Implemented: Binary module renamed to `_plantbox` with a thin `plantbox/__init__.py` wrapper that re-exports the API.
    - Caveats: Dev-source testing copies the built `_plantbox*.so` into `plantbox/` for import; installed wheels will place it correctly via install rules.
  - [x] Add minimal Python package skeleton:
    - Implemented `plantbox/__init__.py` wrapper.
    - `py.typed` deferred.
  - [x] Add CMake install rules for the extension so skbuild can package it correctly (e.g., `install(TARGETS _plantbox ...)`).
  - [x] Populate project metadata in `pyproject.toml` (name, description, license, classifiers, urls, authors).
    - Details: Filled authors, classifiers for CPython 3.9–3.12, URLs (homepage/repo/issues/docs), license file, readme.
  - [x] Decide on versioning strategy:
    - Chosen (for initial wheels): static version `2.1.0` in `pyproject.toml` with `plantbox/_version.py`.
    - Rationale: Our attempt to use SCM-based dynamic versioning via scikit-build-core’s metadata provider and `setuptools_scm` caused build failures in the isolated wheel env (provider requires enabling scikit-build-core experimental plugins). To keep builds green and predictable in Docker without extra flags, we pinned a static version for now.
  - [x] Add a local wheel build helper script (e.g., `scripts/build_wheel.sh`) that builds via `python -m build --wheel` in a clean env/container and prints the artifact path.
  - [x] Add an install-and-smoke helper (e.g., `scripts/test_built_wheel.sh`) that creates a temp venv, installs the wheel, and runs the headless smoke and curated tests from the installed package.
  - [x] Ensure ABI-correct naming is handled by skbuild/pybind11 (remove the temporary aliasing step from the Ubuntu test script once wheels are used for testing).
    - Details: `scripts/run_ubuntu_tests.sh` now builds a wheel and tests the installed package by default; legacy source-tree aliasing removed.
  - [x] Confirm `pip install -e .` works for developer editable installs with scikit-build-core; document any caveats.
    - Details: Verified inside Ubuntu container (venv) with a headless import/simulation; works as expected.
  - [x] Aggregate Ubuntu scripts under a single entry point.
    - Details: Added `scripts/ubuntu/run_all.sh` to run source-tree smoke and wheel smoke sequentially.
  - [x] Add a simple release helper to tag and build artifacts.
    - Details: `scripts/release/release_tag.sh` validates tree, (optionally) runs checks, creates annotated tag (e.g., `v2.1`), optionally pushes, and copies wheels to `release_artifacts/<tag>/ubuntu/`.
  - [x] Ensure wheel artifact versions are correct (no `0.0.0`).
    - Details: Wheels now build and install as `cplantbox-2.1.0-...whl` in the Ubuntu container. We temporarily removed SCM-based dynamic versioning to avoid experimental plugin requirements.

  Blocking tasks discovered (new):
  - [x] Make out-of-source builds find vendored static libs (SUNDIALS/SuiteSparse):
    - Implemented Option A: `src/CMakeLists.txt` now falls back to `${CMAKE_SOURCE_DIR}/src/external/...` when the build-dir copies don’t exist, and the SuiteSparse static link order was corrected. Also build STATIC for wheels via `-DSKBUILD=ON`. Wheel builds now link cleanly in the Ubuntu container.

### Phase 3 — Data layout and import stability

- [x] Package `modelparameter/**` access for installed wheels.
  - Details: Structural subset is installed into the wheel via CMake install rules when building with `SKBUILD`. Wheels can load structural plant/rootsystem XMLs out-of-the-box.
  - Follow-ups: Curate/extend included subsets and keep wheel size reasonable (see below).
- [x] Add a helper `plantbox.data_path()` to get the packaged data root.
  - Details: Helper resolves to package data when installed and falls back to repo tree in dev/source mode.
- [x] Compiled module organization (`_plantbox` + `plantbox/__init__.py`).
  - Details: Completed in Phase 2; wrapper re-exports C++ API and exposes helpers (data path, version).
- [ ] Update a small subset of tutorials/tests to use the helper instead of fragile relative paths (without breaking legacy usage where reasonable).

  Additional Phase 3 tasks (new):
  - [x] Add a tiny unit test ensuring `plantbox.__version__` exists and follows PEP 440 (dev/local builds allowed).
    - Implemented in `test/test_version_and_data.py`.
  - [x] Implement and test `plantbox.data_path()` to work both from source tree and installed wheel.
    - Implemented in `plantbox/__init__.py`; tested in `test/test_version_and_data.py`.
  - [ ] Curate packaged data subsets; consider switching to `tool.scikit-build` packaging config once stable.
    - Current: packaging is done via CMake `install(DIRECTORY ...)` during wheel builds (structural subset). Revisit to include only essential files and/or move to TOML-based packaging once stable.
  - [ ] Document data access pattern migration in tutorials (one or two examples updated first).

### Phase 4 — External dependency strategy for portable wheels

- [ ] Add CMake options to select dependency strategy:
  - `USE_SYSTEM_SUITESPARSE` / `USE_SYSTEM_SUNDIALS` (default OFF for wheels).
  - `BUNDLE_SUITESPARSE` / `BUNDLE_SUNDIALS` to build minimal static libs (AMD/COLAMD/BTF/KLU; SUNDIALS CVODE + required sun*).
- [ ] Validate redistribution licenses for bundled components.
- [ ] Ensure link order and symbols resolve without pulling external BLAS/LAPACK.
- [ ] Verify Linux manylinux compliance (static or vendored non-allowed libs) by auditing built wheels.

### Phase 5 — Multi-platform wheels (local first)

- [ ] Add `cibuildwheel` config to build wheels locally:
  - Linux: manylinux2014 x86_64 and aarch64.
  - macOS: x86_64 and arm64 (or universal2), set `CMAKE_OSX_DEPLOYMENT_TARGET`.
- [ ] Use `auditwheel`/`delocate` to finalize binaries; run our headless smoke tests against produced wheels.
- [ ] Document local build instructions for maintainers.

### Phase 6 — Pybind11 version strategy

- [ ] Stabilize on a known-good `pybind11` (2.11.x) for initial wheel releases.
- [ ] Evaluate/port to pybind11 3.x once wheels are green:
  - Update includes and Python discovery.
  - Re-run local `cibuildwheel` builds and smoke tests.

### Phase 7 — Developer UX and docs

- [ ] Add quickstart docs for devs: build from source, run smoke tests, build wheels locally.
- [ ] Optionally introduce `uv` for fast Python envs (development only), while keeping packaging via `pyproject`.
- [ ] Document test tiers: curated smoke vs. full (including VTK) for capable environments.

---

## Risks & challenges

- SuiteSparse/SUNDIALS static builds across platforms/architectures; wheel sizes and license compliance.
- VTK and GUI dependencies make tests/tutorials non-headless by default; we will keep a headless smoke path and mark GUI tests.
- Packaging large resource trees (`modelparameter`) while keeping installs lean and paths robust.

## Success criteria

- Headless smoke test is green on every change.
- `pip install` wheels work on Linux and macOS for CPython 3.9–3.12 (initially), without extra system deps.
- Tutorials that don’t require GUI run from the installed package with packaged data.

---

### Notes / Changelog

- 2025-08-08: Added Ubuntu amd64 docker test env + runner; curated headless smoke passing; wrote this plan.
- 2025-08-08: Phase 2 wheel builds are green; pinned static version `2.1.0` for initial wheels to avoid scikit-build-core experimental metadata provider; wheel artifacts now named `cplantbox-2.1.0-...whl`.

### Appendix — Versioning note (why static now, how to reintroduce SCM later)

- Why static now
  - scikit-build-core’s metadata “provider” mechanism for dynamic versioning raised: "experimental must be enabled currently to use plugins not provided by scikit-build-core" in the isolated build container.
  - We want deterministic, frictionless local wheel builds without enabling experimental flags or relying on VCS metadata in the build env.
  - Therefore, we set `project.version = "2.1.0"` and keep a simple `plantbox/_version.py` so `plantbox.__version__` is stable.

- How to reintroduce SCM later (two options)
  1) Classic `setuptools_scm` (no scikit-build-core metadata plugin):
     - In `pyproject.toml`: remove the static `project.version`, set `project.dynamic = ["version"]`, add `setuptools_scm[toml]>=8` to `[build-system].requires`, and configure `[tool.setuptools_scm] write_to = "plantbox/_version.py"`.
     - Ensure Git tags are available in the build context (avoid shallow/metadata-less clones). Optionally set `SETUPTOOLS_SCM_PRETEND_VERSION` in release scripts for reproducibility.
  2) scikit-build-core metadata provider (once stable):
     - Use `[tool.scikit-build.metadata] version.provider = "setuptools_scm"` with the appropriate non-experimental configuration (when available); drop the static version.
     - Keep our release helper to create annotated tags; dynamic versions derive from tags.
