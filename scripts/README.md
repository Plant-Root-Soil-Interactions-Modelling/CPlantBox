scripts/

Quick entry points for local builds/tests. Legacy flows are preserved under `old_scripts/` for reference.

- Ubuntu (transitional guardrail; will be superseded by manylinux smoke):
  - Build image: `scripts/ubuntu/build_image.sh`
  - Build+test wheel: `scripts/ubuntu/smoke_wheel.sh`
  - Orchestrator: `scripts/ubuntu/run_all.sh`

- Common helpers (sourced by scripts):
  - `scripts/common/env.sh`: deterministic env setup
  - `scripts/common/logging.sh`: logging helpers
  - `scripts/common/docker.sh`: buildx/image/run helpers

Notes

- Requires Docker Desktop with buildx/QEMU to run Ubuntu amd64 on Apple Silicon.
- Exports `OMP_NUM_THREADS=1` for deterministic results.
