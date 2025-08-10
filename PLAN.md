## CPlantBox Packaging & Wheels — Next Steps Plan (2025)

### Context and current state

Over the last iteration, we modernized packaging and stabilized cross‑platform wheels:

- Packaging/build
  - Adopted `pyproject.toml` with `scikit-build-core`; green wheel builds on Linux manylinux and macOS (local).
  - CMake hygiene: modern Python discovery, safe defaults (Release), portable flags, IPO/LTO disabled for container reliability, corrected SuiteSparse link order.
  - Linux manylinux wheels (x86_64, aarch64) build and repair via `auditwheel`; artifacts smoke‑tested.
- Data and import stability
  - Packaged structural `modelparameter/**` subset; added `plantbox.data_path()`; tests/examples migrated selectively.
- Golden headless test
  - Unified runner `scripts/common/golden_test.py`; standardized render geometry; tolerant similarity (>= 0.9); preserves generated and reference images.
  - Reliable headless execution: manylinux wheel smoke is executed in an Ubuntu headless container (Xvfb + Mesa) for stability; VTK pinned per arch (amd64: 9.2.6, arm64: 9.5.0).
- Scripts reorganization
  - New tree under `scripts/` with platform subfolders and common helpers (`env.sh`, `docker.sh`, `logging.sh`).
  - Manylinux build/test in `scripts/manylinux/{build,test}`; Ubuntu guardrail in `scripts/ubuntu/{build,test}`; macOS `build/` and `test/` in place; legacy scripts archived.

Result: We can build, repair, and run the headless golden test for Linux manylinux wheels (CPython 3.9–3.12 on x86_64 and aarch64), and validate macOS wheels on the host via the same golden. Similarity scores are stable (~0.90–0.95).

---

## Plan — Phases and tasks

We focus on clear, shippable steps. Check off tasks as they complete; keep notes succinct.

### Priorities and execution order

- Priority: macOS developer and wheel flows first (fast local iteration on Apple Silicon host).
- Then: shared hardening (golden polish, size/strip) and Linux stabilization.
- CI: bring up macOS lanes first, then add Linux manylinux lanes.

### Phase 1 — macOS completion (PRIORITY)

Why this phase: Our primary development machine is an Apple Silicon Mac, and macOS builds run fastest locally. Completing the macOS loop first maximizes iteration speed, gives immediate developer feedback, and delivers a high‑quality user experience for the largest contributor base before we generalize elsewhere.

- [ ] `scripts/macos/run_all.sh` orchestrator: editable smoke → curated/golden → build/repair → per‑interpreter smoke.
- [ ] `scripts/macos/dev_editable.sh` (new) using Homebrew deps; fast inner‑loop.
- [ ] `scripts/macos/run_tutorials_headless.sh` with curated headless‑safe examples using `plantbox.data_path()`.
- [ ] Delocation diagnostics and policy
  - `delocate-listdeps --all` preflight + fail‑fast on missing/non‑system deps.
  - RPATH audit on repaired `_plantbox*.so`.
  - Autodetect/enforce `MACOSX_DEPLOYMENT_TARGET` (arm64 default 15.0 unless overridden); detect Brew min‑OS of linked dylibs and warn/fail if mismatched.
- [ ] Optional: universal2 evaluation (or stay per‑arch with clear guidance).
- [ ] Speed‑ups for macOS local builds
  - Reuse a persistent venv for smoke; skip reinstall if unchanged.
  - Ensure Ninja is used; set `CMAKE_OSX_ARCHITECTURES` explicitly; pass `-DCMAKE_BUILD_TYPE=Release`.
  - Gate expensive steps behind flags (e.g., `--rebuild-third-party`).

### Phase 2 — Stabilization and hardening (Linux + shared)

Why this phase: After fast macOS iteration, we need cross‑platform rigor. Finalizing dependency switches ensures maintainers can pick system vs bundled libs safely; size/strip keeps artifacts lean; golden polish reduces flakes across headless environments; script UX and manifests improve operability and traceability. This phase converts “works on my machine” into reproducible, auditable outputs on all Linux targets while keeping shared tooling consistent.

- [ ] Finalize dependency switches in CMake
  - `USE_SYSTEM_SUITESPARSE` / `USE_SYSTEM_SUNDIALS` (default OFF for wheels)
  - `BUNDLE_SUITESPARSE` / `BUNDLE_SUNDIALS` (default ON for wheels)
  - Enforce defaults automatically for `SKBUILD=ON`; document overrides.
- [ ] Wheel size/strip and audits
  - Add strip step and print size summaries per repaired wheel.
  - Always run `auditwheel show` and capture external deps in logs.
- [ ] Golden test polish
  - Preflight diagnostics if golden asset missing; actionable guidance.
  - `CPB_VERBOSE=1` gate for extra logs.
  - Optional: `CPB_HEADLESS_LINE_MODE=1` rendering path to further reduce aliasing.
- [ ] Script UX
  - Add `-h/--help` to key scripts; consistent error messages.
  - Document env knobs in `scripts/README.md` (`CPB_GOLDEN_OS`, `CPB_GOLDEN_THRESH`, `CPB_DOWNSCALE`, `CPB_DIFF_THRESH`).
- [ ] Wheelhouse manifests
  - Emit `index.txt` in each wheelhouse subdir with timestamp, git describe, python tag(s), tool versions.

Notes:

- Keep default golden enabled in smokes; continue using Ubuntu headless container to execute manylinux golden.

### Phase 3 — cibuildwheel + CI (start with macOS)

Why this phase: Once macOS is reliable and shared tooling is hardened, automating builds/tests ensures repeatability and broad coverage. Starting with macOS leverages the work in Phase 1 and brings fast signal in CI; adding Linux next extends coverage to our distribution targets.

- [ ] Add cibuildwheel configuration
  - Linux: manylinux2014 `x86_64 aarch64`, CPython 3.9–3.12.
  - macOS: `x86_64 arm64` (or universal2 if chosen); repair with `delocate`.
  - Repair commands: Linux `auditwheel repair`, macOS `delocate-wheel`.
- [ ] GitHub Actions (or CI of choice)
  - Start with macOS lanes (arm64, x86_64), cache, artifact upload.
  - Add Linux lanes next.
  - Post‑install golden smokes:
    - macOS: run unified golden on the runner host.
    - Linux: run unified golden inside our Ubuntu headless container against the produced manylinux wheel.
  - Tag workflow to publish release artifacts (manual for now).
- [ ] Nightly/cron smoke job (latest main): build one wheel per platform and run golden smokes.

### Phase 4 — Developer UX and docs

Why this phase: Great tooling without guidance slows adoption. Consolidated docs and helper scripts reduce onboarding time, standardize workflows, and improve reliability across contributors.

- [ ] Expand `scripts/README.md` with quickstarts per platform and env knobs.
- [ ] Top‑level `scripts/run.sh` dispatcher for common tasks.
- [ ] Tutorial migration
  - Extend migration to `plantbox.data_path()` where appropriate; retain legacy compat.
  - Curate packaged data to keep wheels lean; consider TOML packaging once stable.
- [ ] Contributor docs: build from source, smokes, wheels, troubleshooting.

### Phase 5 — Versioning and releases

- [ ] Reintroduce SCM versioning (preferred: classic `setuptools_scm` writing `plantbox/_version.py`).
- [ ] Update release helper(s): tag, build matrix locally/CI, verify, publish.
- [ ] Keep CHANGELOG snippet in releases (auto‑generated or manual).

### Phase 6 — Typing and IDE experience (macOS dev loop benefit)

Why this phase: The compiled extension hides symbols from static tooling, degrading IDE assistance. Shipping stubs and richer pybind11 metadata makes the API discoverable, improves code completion and inline docs, and reduces user friction—especially impactful for the macOS dev loop we prioritize.

Problem: IDEs and static analyzers do not see symbols like `Plant` on `plantbox` because the API is defined in a C++ extension. We need to ship stubs and improve introspection.

- [ ] Publish `.pyi` stubs with the wheel
  - Generate initial stubs using `pybind11-stubgen` (CI/build step) or `stubgen`, then curate a minimal public surface (e.g., `Plant`, `SegmentAnalyser`, key functions).
  - Place curated stubs under `plantbox/` (e.g., `plantbox/__init__.pyi` plus module stubs) so editors and type checkers resolve attributes.
  - Ensure stubs are included in wheels via packaging config (CMake install or `tool.scikit-build` include).
- [ ] Improve runtime introspection/docstrings
  - In pybind11 bindings, add `py::arg("name")` and docstrings to functions/constructors; disambiguate overloads with `py::overload_cast` to expose clean signatures.
  - Expose `__all__` on the module; keep thin Python shim `plantbox/__init__.py` re‑exporting the C++ module to centralize docs and helpers.
- [ ] Optional: mark package typed
  - If we add inline type hints in the Python shim, include `plantbox/py.typed` (PEP 561). Not required if we rely solely on `.pyi` stubs.
- [ ] Validation
  - Verify symbol visibility and signatures in VS Code (Pylance) and PyCharm.
  - Add a minimal `mypy` job (non‑blocking) to sanity‑check stubs for basic usage.

---

## Quick commands (current)

- Linux manylinux wheels (build & smoke):
  - Build: `./scripts/manylinux/build/build-wheels.sh x86_64` (or `aarch64`)
  - Smoke: `./scripts/manylinux/test/test-wheel.sh x86_64` (reuses Ubuntu headless golden)
- Ubuntu smoke (guardrail): `./scripts/ubuntu/run_all.sh`
- macOS wheel validation (host): `./scripts/macos/test/test-wheel.sh`

Guardrails:

- Golden tests are enabled by default; similarity threshold 0.9.
- Determinism: export `OMP_NUM_THREADS=1` in scripts; pinned VTK in containers.
