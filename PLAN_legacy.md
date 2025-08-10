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
  - Note: This Ubuntu image is a convenience for fast, repeatable local smoke tests. Our Linux distribution target is manylinux wheels (portable across glibc-based distros). As we add manylinux smoke tests, the Ubuntu guardrail may be retired.

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
- [x] Update a small subset of tutorials/tests to use the helper instead of fragile relative paths (without breaking legacy usage where reasonable).
  - Details: Updated `tutorial/chapter2_structure/example_plant.py` and `tutorial/chapter2_structure/example_plant_anim.py` to use `plantbox.data_path()`; wheel smoke test now loads parameters via `data_path()` as well.

  Additional Phase 3 tasks (new):
  - [x] Add a tiny unit test ensuring `plantbox.__version__` exists and follows PEP 440 (dev/local builds allowed).
    - Implemented in `test/test_version_and_data.py`.
  - [x] Implement and test `plantbox.data_path()` to work both from source tree and installed wheel.
    - Implemented in `plantbox/__init__.py`; tested in `test/test_version_and_data.py`.
  - [ ] Curate packaged data subsets; consider switching to `tool.scikit-build` packaging config once stable. (Moved to Backlog)
    - Current: packaging is done via CMake `install(DIRECTORY ...)` during wheel builds (structural subset). Revisit to include only essential files and/or move to TOML-based packaging once stable.
  - [ ] Document data access pattern migration in tutorials (expand beyond the two examples; keep legacy paths working where feasible). (Moved to Backlog)

### Phase 4 — External dependency strategy for portable wheels

- [ ] Add CMake options to select dependency strategy:
  - `USE_SYSTEM_SUITESPARSE` / `USE_SYSTEM_SUNDIALS` (default OFF for wheels).
  - `BUNDLE_SUITESPARSE` / `BUNDLE_SUNDIALS` to build minimal static libs (AMD/COLAMD/BTF/KLU; SUNDIALS CVODE + required sun*).
- [ ] Validate redistribution licenses for bundled components. (Moved to Backlog)
- [x] Ensure link order and symbols resolve without pulling external BLAS/LAPACK.
  - Details: Corrected SuiteSparse link order earlier and embedded static `CPlantBox` via whole-archive for Linux wheel builds; repaired manylinux wheels import/run without external BLAS/LAPACK.
- [x] Verify Linux manylinux compliance (static or vendored non-allowed libs) by auditing built wheels.
  - Details: Built in `manylinux2014_x86_64` container and ran `auditwheel show/repair`. Produced repaired wheels for CPython 3.9–3.12; in-container import+simulate smoke PASS.

  Readiness: Yes. Linux wheels are built as manylinux-compliant artifacts and run green; structural data packaged; versioning stable.

  Additional Phase 4 tasks (added):
  - [x] Build and repair manylinux wheels (x86_64) for CPython 3.9–3.12.
    - Details: `scripts/wheels/build_manylinux2014.sh` builds with container Pythons, sets CMake Python hints, runs `auditwheel repair`, and smoke-tests each wheel.
  - [x] Add a wheel audit step (Ubuntu for now): run `auditwheel show` on built wheels and capture external deps (e.g., `libstdc++`, `libgomp`).
    - Details: `scripts/wheels/audit_linux.sh` added; ad-hoc audits also run during manylinux flow.
  - [x] If OpenMP is used, confirm `libgomp` handling (vendor via `auditwheel repair` or disable OpenMP for wheels if not needed).
    - Details: No OpenMP linkage detected in current builds; `libgomp` not present in audit outputs. Keep monitoring if OpenMP is introduced.
  - [ ] Strip and size-audit artifacts; keep the extension self-contained (no RPATH; already stripped on install).
  - [x] Document the dependency strategy matrix (bundled vs. system) and defaults in README/PLAN.
    - Details: Added CMake switches and defaults below; wheels default to bundled. See matrix and override examples.

  Focus now:
  - Implement CMake dependency switches (`USE_SYSTEM_*`, `BUNDLE_*`) and wire system-vs-bundled detection.
  - Add strip and size-audit step for repaired wheels.
  - Document the dependency strategy (matrix + defaults) in README/PLAN.

  Postponed:
  - Licensing/attributions (moved to Backlog).
  - Other platforms (manylinux aarch64, macOS, Windows) – Phase 5+.
  - Broader docs/CI polish – Phase 5+/7.

  Dependency strategy matrix and switches
  - Switches (CMake options in `src/CMakeLists.txt`):
    - `USE_SYSTEM_SUITESPARSE` (OFF by default)
    - `USE_SYSTEM_SUNDIALS` (OFF by default)
    - `BUNDLE_SUITESPARSE` (ON by default)
    - `BUNDLE_SUNDIALS` (ON by default)
    - For wheel builds (`SKBUILD=ON`), defaults force bundled unless explicitly overridden.
  - Source build (developer machines):
    - Default: bundled (vendored static libs). You can opt into system libs:
      Example: `cmake -S . -B build -DUSE_SYSTEM_SUITESPARSE=ON -DUSE_SYSTEM_SUNDIALS=ON -DBUNDLE_SUITESPARSE=OFF -DBUNDLE_SUNDIALS=OFF`
  - Wheel build (packaging):
    - Default: bundled static libs (current tested path). Manylinux wheels are repaired and self-contained.
  - Notes:
    - SuiteSparse static link order is enforced (KLU, AMD, COLAMD, BTF, suitesparseconfig).
    - No OpenMP (`libgomp`) detected in current builds.

### Phase 4.4 — Golden headless tests stabilization

Focus: make the visual regression (golden) test deterministic and portable across headless environments, and ensure it runs reliably for manylinux wheels (using an Ubuntu headless container as execution host when needed).

- [x] Tolerant image similarity metric to handle antialiasing differences
  - Details: `test/tools_image.py::compare_images_png` now supports:
    - Ignoring white-on-white background pixels (as before)
    - Optional 2× block downscaling (`CPB_DOWNSCALE`, default 1)
    - Small per-channel threshold (`CPB_DIFF_THRESH`, default 10)
    - Similarity computed as `1 - MAD` over significant differences
  - Threshold: tests require `similarity >= 0.9` (updated in `test/test_golden_headless.py` and enforced in `scripts/common/golden_test.py`).

- [x] Standardize render geometry to match golden
  - Details: width=1000, height=int(1000*2.2), zoom=3.0 in both `test/test_golden_headless.py` and `scripts/common/golden_test.py`.

- [x] Keep both generated and reference images for inspection
  - Details: Golden smoke copies `test/golden/<os>/example_plant_headless.png` to `test/golden/example_plant_headless.png` and writes `test/golden/generated_by_test.png`; both persist after the run.

- [x] Reliable headless execution environment
  - Details: Golden run executes inside the Ubuntu headless container under `xvfb-run` with `LIBGL_ALWAYS_SOFTWARE=1` and VTK pinned (`9.2.6` on amd64, `9.5.0` on arm64) to avoid blank frames.
  - Scripts: `scripts/ubuntu/test/test-wheel.sh` (Ubuntu wheels) and `scripts/manylinux/test/test-wheel.sh` (manylinux wheels) install the wheel in a container venv and execute `scripts/common/golden_test.py`.
  - Note: manylinux in-container headless rendering remains flaky; our manylinux smoke reuses the Ubuntu headless runner to ensure stable offscreen rendering.

- [ ] Optional: line rendering mode to further reduce aliasing
  - Task: Implement `CPB_HEADLESS_LINE_MODE=1` in `test/tools_image.render_headless_png` (render polyline actors instead of tubes) and enable it in `smoke_on_ubuntu.sh` for golden runs.

- [ ] Golden assets preflight and diagnostics
  - Task: Add explicit preflight in `golden_smoke.py` to validate golden presence and emit actionable guidance if missing (partial checks exist; improve messages and help hints).

- [ ] Script verbosity control
  - Task: Gate extra diagnostics behind `CPB_VERBOSE=1` to keep default logs concise, while retaining a detailed mode for debugging.

- [ ] Document env knobs in scripts/README.md
  - Task: Document `CPB_DOWNSCALE`, `CPB_DIFF_THRESH`, `CPB_HEADLESS_LINE_MODE`, and how to tweak thresholds when investigating golden divergences.

### Phase 4.5 — Tidying up scripts folder

Streamline build/test entry points while keeping the latest Ubuntu golden-flow intact (transitional). We moved existing scripts to `old_scripts/` and constructed a clean `scripts/` tree. Legacy remains for bookkeeping. Linux distribution targets manylinux; Ubuntu-specific runners are for local smoke and also serve as the golden host for manylinux wheels.

- [x] Canonicalize Linux guardrail around manylinux (Ubuntu is transitional):
  - Transitional: `old_scripts/run_ubuntu_tests.sh` remains the reference headless smoke with golden test until a manylinux-based smoke is in place. (kept, unchanged)
  - New: Added manylinux wheel build and smoke flows and a dedicated Ubuntu-based golden smoke for manylinux wheels:
    - Build/repair: `scripts/manylinux/build/build-wheels.sh` (per-arch builder image; `auditwheel repair`)
    - Golden smoke runner: `scripts/manylinux/test/test-wheel.sh` selects a matching wheel and runs golden via the Ubuntu headless container.
  - Note: Ubuntu golden-runner is currently canonical for Linux wheels; manylinux containers alone are not reliable for headless VTK.
  - [ ] Add a short usage note in `scripts/README.md` explaining the transition from Ubuntu to manylinux guardrails.

- [x] Introduce common helpers for reuse:
  - Implemented `scripts/common/env.sh` (deterministic env; sets `OMP_NUM_THREADS=1`), `scripts/common/logging.sh` (uniform logs), and `scripts/common/docker.sh` (buildx helpers).
  - New scripts source these helpers. Legacy scripts are preserved and not modified.

- [x] Ubuntu scripts (transitional; legacy preserved in `old_scripts/`):
  - `scripts/ubuntu/build_image.sh` builds `docker/Dockerfile.ubuntu-test-env` (linux/amd64 via buildx).
  - `scripts/ubuntu/smoke_wheel.sh` builds a wheel, installs in a temp venv, runs a minimal simulation and the golden test.
  - `scripts/ubuntu/run_all.sh` orchestrates the two.
  - Marked as convenience runners; will be deprecated once manylinux smoke reaches parity.

- [ ] macOS local developer flow (new implementation; legacy preserved in `old_scripts/`):
  - Implement `scripts/macos/run_all.sh` to orchestrate: editable install smoke → curated pytest subset/golden headless test → build wheels → install & smoke each repaired wheel. (in progress)
  - Implement `scripts/macos/build/build-wheel.sh` to build wheels with safe defaults and clear artifact output. (done)
  - Implement `scripts/macos/test/test-wheel.sh` to venv-install a macOS wheel, install VTK, and run the unified golden. (done)
  - Implement `scripts/macos/run_tutorials_headless.sh` to execute a curated list of headless-safe tutorials that use `plantbox.data_path()` and run from arbitrary CWD. (todo)
  - Implement `scripts/macos/dev_editable.sh` anew, based on the legacy logic but adapted to the new helper structure. (todo)

- [x] Wheels scripts (new implementation; legacy preserved in `old_scripts/`):
  - Implemented `scripts/manylinux/build/build-wheels.sh` and per-arch Dockerfiles to build/repair manylinux wheels.
  - Implemented `scripts/manylinux/test/test-wheel.sh` to smoke-test repaired manylinux wheels using the Ubuntu headless runner.
  - Legacy `scripts/wheels/*` moved to `old_scripts/` and replaced.

- [x] Golden image test consistency:
  - Unified golden test script at `scripts/common/golden_test.py` used by all runners. It copies and compares the correct golden (`test/golden/<os>/example_plant_headless.png` → `test/golden/example_plant_headless.png`) and keeps both images for inspection.
  - [ ] Add/strengthen preflight to verify golden assets exist and print actionable guidance if missing.

- [x] Script UX and docs:
  - Added `scripts/README.md` with overview and usage for new entry points; will expand with env knobs and guardrail notes.
  - Standardized `set -euo pipefail` and deterministic behavior (`OMP_NUM_THREADS=1`) in new scripts; `PYTHONPATH` handling fixed in golden smoke to avoid source-vs-wheel import clashes.
- [ ] Add `-h/--help` flags across new scripts and expand README with env tuning knobs (document `CPB_GOLDEN_OS`, `CPB_GOLDEN_THRESH`, `CPB_DOWNSCALE`, `CPB_DIFF_THRESH`).

- [ ] Hygiene and safety:
  - Run `shellcheck` locally on new/changed scripts and fix high-signal warnings.
  - Ensure `git submodule update --init --recursive` is called where needed in container workflows.
  - Keep old scripts for reference by adding the `old_` prefix instead of deleting. (done earlier)

Migration to the “ideal” scripts layout (step-by-step)

- [x] Step 1 — Preserve current state, then create a clean slate:
  - Rename the current `scripts/` directory to `old_scripts/` (done).
  - Create a fresh `scripts/` directory with the ideal layout skeleton and placeholders (done; now populated with working flows).
  - Add `scripts/README.md` with a quickstart and a note that legacy flows remain in `old_scripts/` (done; expand with env knobs next).

- [x] Step 2 — Create the ideal layout skeleton (empty or minimal scripts that print help):

  Target structure:

  ```
  scripts/
    README.md
    run.sh
    common/
      env.sh
      docker.sh
      python.sh
      logging.sh
    ubuntu/
      build_image.sh
      smoke_wheel.sh
      run_all.sh
    macos/
      dev_editable.sh
      build_wheels.sh
      run_tutorials_headless.sh
      run_all.sh
    wheels/
      linux/
        manylinux_x86_64.sh
        manylinux_aarch64.sh
        audit.sh
      macos/
        build.sh
        delocate_audit.sh
      build_all_local.sh
    release/
      release_tag.sh
  ```

- [x] Step 3 — Port the Ubuntu guardrail into the new layout:
  - `ubuntu/build_image.sh`, `ubuntu/smoke_wheel.sh`, and `ubuntu/run_all.sh` implemented; golden is executed in-container. Parity validated via passing smokes.

- [ ] Step 4 — Port macOS developer and wheel flows:
  - Implement `macos/dev_editable.sh` using logic from `old_scripts/macos/dev_editable.sh` (keep editable flow minimal and reliable).
  - Implement `macos/build_wheels.sh` as a wrapper that calls the current macOS wheel build logic (initially from `old_scripts/wheels/build_macos.sh`), then adds repair and smoke.
  - Add `macos/run_tutorials_headless.sh` with a curated, headless-safe set of examples that use `plantbox.data_path()` and run from arbitrary CWD.
  - Implement `macos/run_all.sh` to orchestrate editable smoke → curated/golden → wheel build/repair → wheel smoke.

- [x] Step 5 — Port Linux wheels (manylinux) flows:
  - Implemented `scripts/manylinux/build/build-wheels.sh` (per-arch) and smoke via `scripts/manylinux/test/test-wheel.sh` using Ubuntu headless runner. Legacy kept under `old_scripts/` for reference.

- [ ] Step 6 — Common helpers and consistency:
  - Implement `common/env.sh` (deterministic env, `OMP_NUM_THREADS=1`), `common/docker.sh` (buildx builder ensure/build/run helpers), `common/python.sh` (temp venv helpers, pip install, wheel install+smoke), and `common/logging.sh` (uniform logs, error handling).
  - Refactor the new scripts to source these helpers; avoid duplication.

- [ ] Step 7 — Wheelhouse layout enforcement and manifests:
  - Ensure all new wheel scripts write to a normalized structure:
    - `wheelhouse/linux/manylinux2014_x86_64/` and `wheelhouse/linux/manylinux2014_aarch64/`.
    - `wheelhouse/macos/arm64/` and `wheelhouse/macos/x86_64/`.
  - Emit an `index.txt` (manifest) per subdir with build metadata (python tag, deploy target, timestamp, git describe, tool versions). (todo)

- [ ] Step 8 — Documentation and references:
  - Update PLAN and README to reference the new `scripts/` entry points (`run.sh`, `ubuntu/run_all.sh`, `macos/run_all.sh`, `wheels/build_all_local.sh`).
  - Clearly state that `old_scripts/` is retained for bookkeeping and should not be modified.

### Phase 5 — Multi-platform wheels (local first)

- [ ] Add `cibuildwheel` config to build wheels locally:
  - Linux: manylinux2014 x86_64 and aarch64.
  - macOS: x86_64 and arm64 (or universal2), set `CMAKE_OSX_DEPLOYMENT_TARGET`.
- [x] Use `auditwheel`/`delocate` to finalize binaries; run our headless smoke tests against produced wheels.
  - Details: manylinux wheels are repaired via `auditwheel` and smoked using the Ubuntu headless runner; macOS wheels are repaired via `delocate` and smoke-tested per interpreter on host using the unified golden.
- [x] Document local build instructions for maintainers (see "Continuous local testing for developers" below).

  Additional Phase 5 tasks (added):
- [x] Linux aarch64: add manylinux aarch64 flow mirroring x86_64 (platform `linux/arm64`), run `auditwheel repair`, and smoke-test wheels using the Ubuntu headless container.
  - Status: Implemented and validated; smokes pass with similarity ≥ 0.9.
  - [x] macOS builder script: add `scripts/wheels/build_macos.sh` to build wheels for CPython 3.9–3.12 on host, with:
    - `CMAKE_OSX_DEPLOYMENT_TARGET` defaulting to 11.0 for arm64 (adjustable via env).
    - `delocate-wheel` to bundle non-system deps.
    - Headless smoke after install for each wheel using packaged data (`plantbox.data_path()`).
    - Status: Build succeeds, but `delocate-wheel` currently fails due to Homebrew deps built with min macOS 14–15, while our default target is 11.0. Workaround: export `MACOSX_DEPLOYMENT_TARGET=15.0` when running the script so the wheel tag matches the bundled dylibs; alternatively rebuild/bundle compatible libs.
    - Note (2025-08-08): Editable install on macOS arm64 verified via `scripts/macos/dev_editable.sh` (Homebrew SUNDIALS/SuiteSparse). Import + minimal simulate passes. Some tutorials still rely on relative paths; they run fine when executed from their directory. Migration to `plantbox.data_path()` will improve this (tracked in Phase 3 / Backlog).
  - [ ] `cibuildwheel` config in `pyproject.toml` (or `cibuildwheel.toml`) with per-platform repair commands:
    - Linux: `CIBW_REPAIR_WHEEL_COMMAND_LINUX="auditwheel repair -w {dest_dir} {wheel}"`.
    - macOS: `CIBW_REPAIR_WHEEL_COMMAND_MACOS="delocate-listdeps {wheel}; delocate-wheel -w {dest_dir} {wheel}"`.
    - Set `CIBW_BUILD` to cp39–cp312; set `CIBW_ARCHS_MACOS="x86_64 arm64"` (or universal2 strategy if feasible).
  - [x] macOS toolchain flags: ensure no OpenMP is required by default; if ever enabled, bundle `libomp` via `delocate`.
    - Details: No OpenMP detected in current builds; keep monitoring and add `libomp` vendoring if enabled later.
  - [ ] RPath policy (macOS): verify `_plantbox` has no unexpected `@rpath` outside system libs after `delocate`.
  - [ ] Size/strip audit: print repaired wheel sizes and extension `.so`/`.dylib` sizes; set a soft budget and track regressions.
  - [ ] Top-level helper: `scripts/wheels/build_all_local.sh` that orchestrates Linux (manylinux x86_64/aarch64) and macOS builds and runs smoke tests.

  New Phase 5 tasks (added 2025-08-08):
  - [ ] macOS: autodetect/enforce a safe `MACOSX_DEPLOYMENT_TARGET` in `scripts/wheels/build_macos.sh`.
    - Details: On arm64 default to 15.0 unless overridden; probe Homebrew-installed SUNDIALS/SuiteSparse with `otool -l` to determine the maximum minimum OS version among linked dylibs, and set the target accordingly. Fail fast with a clear message if mismatch persists.
  - [ ] macOS: add a preflight delocation diagnostics step.
    - Details: Run `delocate-listdeps --all` on the built wheel and print any missing/non-system libs, including their compatibility and minimum OS versions. Gate the build to prevent publishing unrepaired wheels.
  - [ ] macOS: RPATH and dependency audit after repair.
    - Details: Use `otool -l` and `otool -L` on the repaired `_plantbox*.so` to ensure only system or bundled libs remain and no stray `@rpath` entries are present (complements the existing RPath policy item).
  - [ ] macOS: optional universal2 strategy assessment.
    - Details: Evaluate universal2 vs per-arch wheels. If universal2 is chosen, ensure `CMAKE_OSX_ARCHITECTURES="arm64;x86_64"` builds cleanly and delocation collects both slices.
  - [ ] Wheel smoke: add a minimal tutorial-based smoke that runs from an arbitrary working directory.
    - Details: Execute a headless example that loads parameters via `plantbox.data_path()` to validate packaged data access and avoid reliance on repo-relative paths.

#### Continuous local testing for developers

- Quick Ubuntu guardrail (headless smoke, source + wheel):
  - Requires Docker Desktop. Run: `./scripts/ubuntu/run_all.sh`
  - What it does: builds wheel in a clean Ubuntu image, installs it, and runs curated smoke tests.

- Linux manylinux x86_64 wheels (build, repair, smoke):
  - Requires Docker Desktop (QEMU emulation on Apple Silicon is handled automatically).
  - Run: `./scripts/wheels/build_manylinux2014.sh`
  - Artifacts: `wheelhouse/` in repo root. Script prints audit info and runs a minimal simulation.

- Linux manylinux aarch64 wheels (build, repair, smoke):
  - Requires Docker Desktop with `--platform linux/arm64` support (QEMU). Slower on Apple Silicon.
  - Run: `./scripts/wheels/build_manylinux2014_aarch64.sh`

- macOS wheels (arm64 or x86_64 host):
  - Prereqs: Homebrew `sundials` and `suite-sparse`, and Homebrew Pythons 3.9–3.12.
  - Optional: set `MACOSX_DEPLOYMENT_TARGET` (defaults to 11.0). On Intel, consider 10.15 or 11.0.
  - Run: `./scripts/wheels/build_macos.sh`
  - The script repairs with `delocate` and smoke-tests each built wheel using packaged data paths.

- Fast local dev on macOS (editable install using system deps):
  - Run: `./scripts/macos/dev_editable.sh`
  - Good for iterating on C++/Python code; uses Homebrew SUNDIALS/SuiteSparse.

- General tips:
  - Export `OMP_NUM_THREADS=1` for deterministic timings in tests.
  - Start from a clean tree between builds: remove `dist/` and `_skbuild/` (scripts already do this).
  - Inspect wheels: Linux `auditwheel show dist/*.whl`; macOS `delocate-listdeps dist/*.whl`.
  - Size checks: consider `du -h wheelhouse/*` and `otool -L` (macOS) / `ldd` (Linux container) on the extension to verify dependencies.
  - When loading parameters in tests/examples, prefer `plantbox.data_path()` to avoid brittle relative paths.

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

### Backlog (non-blocking)

- Curate packaged data content for wheels; consider TOML-based packaging once stable; ensure wheel size remains reasonable.
- Expand tutorial/doc migration to use `plantbox.data_path()` beyond the two updated examples; keep legacy paths working where feasible.
- Optional developer UX: add macOS aggregator `scripts/macos/run_all.sh` and a top-level cross-platform dispatcher.
- Optional typing: add `py.typed` and improve type hints in the Python shim as we stabilize APIs.
- Reintroduce SCM-based versioning once scikit-build-core provider is non-experimental or use classic `setuptools_scm` flow; adjust release scripts accordingly.

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
