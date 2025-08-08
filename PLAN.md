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
- [ ] Avoid forcing output directories for libs/modules unless required by wheel build tooling.
- Notes: Audit pending; no behavior changes yet detected beyond standard pybind11 defaults.
- [ ] Introduce pytest markers to tag VTK-dependent tests; skip automatically in headless environments.
- [ ] Add `ruff` (conservative config: E,F,I) and `pre-commit` for python files; document opt-in for devs.
- [ ] Add `clang-format` config (do not enforce globally yet; document usage for contributors).
- [x] Disable LTO/IPO to avoid toolchain issues in container builds.
  - Details: Set `CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF` in `src/CMakeLists.txt` to resolve `lto-wrapper`/jobs server errors seen under emulation.
  - Caveats: Minor potential link-time optimization loss; acceptable for reliability in cross-arch container builds.

  New small hygiene tasks (added):
  - [ ] Set an explicit minimum CMake version compatible with modern `FindPython` (e.g., `cmake_minimum_required(VERSION 3.18)` or higher) and document it.
  - [ ] Add a top-level `.dockerignore` to speed Docker builds and avoid copying large artifacts (e.g., `build/`, `__pycache__/`, `.git/`, `tutorial/results/`).
  - [ ] Make test runs deterministic by default in scripts (export `OMP_NUM_THREADS=1`) and document this for local dev.
  - [ ] Ensure vendored `external/pybind11` is present and at a known-good tag; add a short note in docs and a friendly check in the build script if missing.
  - [ ] Add a simple CTest hook that runs the curated pytest subset (`ctest -R smoke`) to make `cmake --build . && ctest` work out-of-the-box.
  - [ ] Ensure tests do not leave artifacts in the repo; add cleanup and ignore rules.
    - Details: Route test-generated files (e.g., `human.xml`, `leaf.xml`, `organ.xml`, `organism.rsml`, `root.xml`, `rs_parameters.xml`, `seed.xml`, `stem.xml`, `test_rootsystem.rsml`) into per-test temp dirs via `tmp_path` and/or chdir fixtures; add autouse teardown to remove legacy outputs; add `.gitignore` patterns as a fallback.
    - Caveats: Some tutorials intentionally write outputs; keep those under `tutorial/results/` (already ignored) and out of repo root.

### Phase 2 — Introduce modern packaging (pyproject + scikit-build-core)

- [ ] Add `pyproject.toml` using `scikit-build-core` as the build backend.
- [ ] Configure extension build via CMake without setup.py; verify `pip wheel .` works in the Ubuntu test container.
- [ ] Keep runtime python deps minimal (e.g., `numpy`), avoid optional heavy ones (`vtk`) in core install.
- [ ] Prove smoke from fresh venv in-container using the built wheel.

### Phase 3 — Data layout and import stability

- [ ] Package `modelparameter/**` as package data; ensure installed wheels can access it.
- [ ] Add a helper (e.g., `plantbox.data_path()`) to get the packaged data root.
- [ ] (Optional but recommended) Rename compiled module to `_plantbox` and add a thin `plantbox/__init__.py` that imports it and exposes helpers (data path, version), de-risking future changes.
- [ ] Update a small subset of tutorials/tests to use the helper instead of fragile relative paths (without breaking legacy usage where reasonable).

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
