# Sensitivity analysis of the hyphal model — Phase 0 + Morris

Framework for a noise test, `dx` convergence and Morris screening on SLURM, for the
Petri dish setup of `petridishsimauto.py` (Oyarte Galvez et al. 2025) with
`tomatoparameters/TwoHyphaePlusBAS.xml` as the nominal model. Fully tested end to end with a toy model
(`dummy_model.py`).

## What you need to adapt

| File | What |
|---|---|
| `config.py` | **factor ranges** (`distTH` still a ±25 % placeholder), possibly `HOURS_HYPHAE`, `SNAPSHOT_HOURS` |
| `slurm/*.sbatch` | activate `~/cpbenv` (additionally needs SALib, scikit-learn, networkx), partition `cpu256`; possibly time/memory |

`model_hook.py` calls `AMFSimulationSetup.makeSimulation` with the parameter file of the run
and stores the quantities per ring and hour as `petridishsimauto.py` does (`raw.npz`).
Changes to the nominal model belong in `TwoHyphaePlusBAS.xml`, not here.

## Everything at once (`slurm/submit_all.sh`)

```bash
cd ~ && bash CPlantBox/experimental/myc/sa_setup/slurm/submit_all.sh            # everything
bash CPlantBox/experimental/myc/sa_setup/slurm/submit_all.sh --no-morris        # look at Phase 0 first
bash CPlantBox/experimental/myc/sa_setup/slurm/submit_all.sh --help             # all options
```

Order with SLURM dependencies: smoke test (3 runs), then the gate (`smoke_test.py --report --gate`).
If the nominal run is not clean, SLURM cancels everything that follows. After that, Phase 0 → Morris,
the 120 h test and the calibration run in parallel. The designs are each created by a small job (`stage.sbatch`),
which then submits the runs (needs `sbatch` from within jobs). Nothing has to be copied, the sbatch files
come directly from the repo.

## Workflow (individual steps)

Like the Petri dish runs, the jobs are submitted from the **home directory**, and
**every Python script runs as a job** (not on the login node). Once:

```bash
cd ~
cp CPlantBox/experimental/myc/sa_setup/slurm/*.sbatch CPlantBox/experimental/myc/sa_setup/slurm/submit_smoke.sh .
cp CPlantBox/experimental/myc/sa_setup/slurm/submit.sh submit_sa.sh
source cpbenv/bin/activate && pip install SALib scikit-learn networkx
mkdir -p slurmOut slurmErr
```

Everything from `~`. Individual Python steps with `sbatch sa_run.sbatch <script> [arguments]`,
their output is in `slurmOut/sa_run-<job>.out`:

```bash
cd ~
sbatch sa_run.sbatch check_config.py         # ranges, nominal values, number of runs, CPU hours

# ---- smoke test: 3 real runs (nominal, all factors low, all high), then a report ----
bash submit_smoke.sh                         # or: bash submit_smoke.sh --quick (12 h after t*)
# report: slurmOut/sa_run-<job>.out -> all OK? runtime and memory -> sim.sbatch

# ---- Phase 0: noise, dx convergence, sample of the factor space (100 runs) ----
sbatch sa_run.sbatch design.py phase0        # wait until the job has finished
bash submit_sa.sh results/design_phase0.csv  # design path relative to sa_setup
# collection and analyze_phase0.py run automatically at the end: slurmOut/sa_collect-<job>.out

# ---- test: nominal point over 120 h (10 runs), are the wave fits then defined often enough? ----
sbatch sa_run.sbatch design.py long          # wait until the job has finished
SIMTIME=12:00:00 bash submit_sa.sh results/design_long.csv
# result: slurmOut/sa_collect-<job>.out (analyze_long.py), windows 60, 90 and 120 h

# ---- Phase 1: Morris (r × (k+1) × R runs) ----
sbatch sa_run.sbatch design.py morris        # wait until the job has finished
bash submit_sa.sh results/design_morris.csv
# analyze_morris.py runs automatically at the end, plots in sa_setup/results/plots/
```

To repeat a step by hand, e.g. after changes to the analysis:
`sbatch sa_run.sbatch analyze_morris.py`.

If CPlantBox lives elsewhere, set `SA_DIR=/path/to/sa_setup` before `sbatch` or `bash submit_sa.sh`.
After changes to the sbatch files in the repo, renew the copies in the home directory.

Testing locally without a cluster and without the real model:

```bash
export SA_DUMMY=1
python3 design.py phase0
python3 run_one.py results/design_phase0.csv 0-4
python3 evaluate.py results/design_phase0.csv 0-4
```

**Aborted runs** (timeout, node failure): just call `submit_sa.sh` again,
finished runs are skipped. `results/status_*.csv` shows what is missing. Runs without t*
after `MAX_WAIT_DAYS` count as failed (error message in `runs/<id>/meta.json`).

**New outputs / other times:** change `evaluate.py` or `config.SNAPSHOT_HOURS`, then
`cd ~ && bash submit_sa.sh <design> --eval-only`. No re-simulation needed, `raw.npz` contains
all quantities per ring and hour. Network quantities (`AMFNetwork.py`) need the
geometry: set `SAVE_VTP = True` for that.

## Factors (13)

| Group | Factors | Range |
|---|---|---|
| Runner (subType 1) | `R_v`, `R_hlt`, `R_lmax`, `R_ln`, `R_theta` | v 0.432–0.72 cm/d (180–300 µm/h), hlt 4–10 d, lmax 1.5–5 cm, ln 0.1–0.6 cm (log), θ 50–80° |
| BAS (subType 2) | `B_v`, `B_b`, `B_hlt`, `B_theta` | v as runner, b 0.55–1.35 /d, hlt 5–9 d (nominal 7 d), θ 50–80° |
| shared | `distTH`, `ana`, `v_cv` | ±25 %, 0.5–1, CV 0.15–0.40 |
| control | `R_b_ctrl` | b of the runners, no effect with b_prob = 0 |

The sources of the ranges are given per factor in `config.py` (`source`).

Fixed: successor split (0.75 BAS / 0.25 runner), la, lb, BAS lmax (≈ 3 × dx),
radius a, tropism, dx, b_prob.

## Outputs

| Kind | Outputs | Time |
|---|---|---|
| quantities | `total_length` [cm], `active_tips`, `bas_share`, `anastomoses` (since t*), `front_position` (tip-weighted radius, mm) | 10, 30, 60 h after t* |
| wave | `front_velocity` [mm/d], `K1` [µm mm⁻²], `K2` [mm⁻²], `wave_width` [h], `encounter_C`, `beta` [µm/h] (a = β ρ n, target ≈ 23) | fits over 0–60 h, stored at 60 h |

The wave quantities come from the fits of `AMFDensities` with the quality criteria of
`petridishstats.py` (R², saturation, at least 3 rings). If the fits fail, the value is missing,
and the run counts like a failed one for this output. The number of rings behind
the fits is stored as diagnostics (`n_wave_rings`, `n_K1_rings`, `n_tip_rings`) in `outputs.csv`,
as is `wait_days` (time until t*). Diagnostics do not enter the SA.

## Design decisions

- **The `value`s are varied.** With `DEV_MODE = "relative"`, `dev` is scaled
  proportionally (constant CV). Exception `v_cv`: sets dev of v to CV × v, for both subtypes.
- **Shared factors** (`distTH`, `ana`, `v_cv`) are written to runner and BAS.
- **b_prob as structure, not as a factor.** b is only read during tip splitting
  (`Hyphae.cpp`), so `R_b_ctrl` is a **control factor**: it must give μ* = 0.
- **No abort on a stalling front** (`MIN_GROWING_TIPS = None`): in the SA a
  slow front is a result, not an error. Runs are only aborted without t*.
- **`dx` only in Phase 0** (convergence test), never in the sensitivity analysis.
- **Shared seeds along a trajectory** (common random numbers): elementary
  effects are differences, equal seeds remove a large part of the noise.
- **Noise level from the replicates:** every replicate is its own consistent
  trajectory. The spread of the EEs between replicates gives, per factor, the μ*
  that pure seed noise would produce. Factors with μ* ≤ 2 × this level count as
  "within noise" (column `in_noise`, black tick in the ranking plot).
- **μ* in normalised units:** output change per traversal of the factor range,
  comparable between factors. `mu_star_rel` additionally normalises per output.
- **Morris via SALib** (`pip install SALib`): design with `SALib.sample.morris`
  (trajectory selection by maximum spread, `local_optimization`), analysis with
  `SALib.analyze.morris` on the replicate mean. Factors are sampled on `[0,1]`
  and converted in `design.to_physical` (also log-uniform). Only the noise level
  is computed here, because SALib does not return the individual EEs.
- **Failed runs:** a replicate of a trajectory only counts if all k+1
  points exist; a trajectory counts if at least one replicate is complete.

## What Phase 0 decides

- `phase0_noise.csv`: CV per output. With a high CV (> ~0.3), rather increase `MORRIS_REPLICATES`.
  Column `frac_zero` shows sparsely populated outputs.
- `phase0_dx.csv`: if the nominal `dx` deviates significantly (> 2 %) from the finest, `dx`
  is too coarse. Check especially for anastomoses (distTH = dx) and BAS (lmax ≈ 3 × dx).
  A small dx costs a lot of runtime: increase `--time` in `sim.sbatch` for it.
- `phase0_dt.csv`: the same for the time step (1 h, 0.5 h, 0.25 h; `config.PHASE0_DT_STEPS_PER_DAY`). The 1 h comes
  from the image interval of the experiment, which only fixes the evaluation interval (stays hourly). Anastomosis is
  checked once per step at the tip, `ana` and the BAS tip splitting act per step: check `beta`, `anastomoses`,
  `encounter_C` and `bas_share`. Small dt costs runtime roughly ∝ steps: increase `--time` in `sim.sbatch`.
- `phase0_status.csv`: how many of the 20 random points in the factor space reach t*, and how
  long they run. Many aborts: ranges too wide, otherwise whole trajectories go missing in Morris.
- Runtime and memory: adjust `--time`/`--mem` in the sbatch files.

## Results

- `results/morris_indices.csv`: μ*, CI, μ, σ, noise level per `time × region × output × factor`
- `results/plots/morris_ranking_t*.png`: μ* with bootstrap CI (SALib: μ* ± half-width) per output
- `results/plots/morris_musigma_t*.png`: μ* vs. σ (σ ≳ μ*: nonlinear/interactions)
- `results/plots/morris_overview.png`: relative sensitivity of all outputs and times

`beispiel_dummy/` contains plots from a test run with the toy model.
They only show the format, their content means nothing.

## Calibration (`calibrate.py`)

Root finding `F(x) = 0` with `F_i = log(y_i / y*_i) / log(1 + tol_i)`. A target counts as
reached when |F_i| ≤ 1. Solved with `scipy.optimize.least_squares` (TRF, box bounds)
in the normalised parameter space. The Jacobian comes from large forward differences (≈ 0.15),
all difference points run at the same time. The search stops as soon as all targets are within tolerance.

16 parameters: separate for runner/BAS as in the SA (v, hlt, theta, b, ln, lmax), plus separate
tropismN and the radius a, and the BAS share of the runner successors (`successor percentage`).
Shared: ana, distTH. lb of the runners is coupled to ln, la stays fixed (template XML).
Current targets: β 23 µm/h, K1 1400 µm mm⁻², BAS share 0.3, K2 0.15 mm⁻² (runner tips), branching rate
4 %/h (new runner tips per growing runner tip and hour, `branch_new_runner` in `raw.npz`),
wave speed 280 µm/h (all Oyarte Galvez et al. 2025) and relative root efficiency
≈ 0.5: (E_r − E_r,MST)/(E_r,DT − E_r,MST) between tips/BAS attachment sites and root nodes (≤ 5 mm
from the barrier), at the end of the simulation (`AMFNetwork.relativeRootEfficiency`, can be switched off via
`config.ROOT_EFFICIENCY`). As in the experiment, a planar graph seen from above is evaluated
(`AMFNetwork.getPlanarNetwork`): all hyphal segments on the fungal side projected onto the x-y plane, every
crossing of two segments becomes a node (degree 4 in the image skeleton).
Because 7 targets do not determine 16 parameters uniquely, every parameter has a **prior**, an
additional residual z = (x − centre)/σ (log for log parameters). Literature priors (`CalParam.prior`, range
= ±2σ) exist for R_v, B_v, R_hlt, R_ln, R_theta, B_theta, B_hlt, B_a (BAS diameter 2–10 µm, Bago 1998) and R_a
(radius 50 µm, Schnepf 2017, prior 25–100 µm; lower bound of the search range 10 µm, Bago 1998). All others have a weak
prior around the nominal value (σ = `WEAK_PRIOR_SIGMA` × search range). Parameters that no target acts on
thus stay at the nominal value. The report flags literature priors from which the solution deviates by more than 2σ.
All evaluations use the same seeds (CRN), the outputs are pooled over the seeds, the
standard error comes from a bootstrap. Parameters, search ranges, targets and tolerances
are at the top of `calibrate.py` (`PARAMS`, `TARGETS`).

```bash
cd ~ && cp CPlantBox/experimental/myc/sa_setup/slurm/calibrate.sbatch .
sbatch sa_run.sbatch calibrate.py --check    # ranges, start point, number of evaluations
sbatch calibrate.sbatch                      # 128 cores, 16 seeds per evaluation
```

Results in `calibration/`: **`report.txt`** (rewritten after every evaluation, also after an abort:
status, best point, targets, failed runs with error messages, undefined outputs, noise,
parameters without influence, Jacobian rank), `evals.csv` (every evaluation with values, SE, F),
`jacobian.csv` (sensitivities dF/du, rewritten after every Jacobian),
`best_params.xml`. Restartable: submit the same command again.
Local test: `SA_DUMMY=1 python3 calibrate.py --seeds 8`.

## Literature

- Morris (1991), *Technometrics* 33:161–174
- Campolongo, Cariboni, Saltelli (2007), *Environ. Model. Softw.* 22:1509–1518 (μ*, trajectory selection)
- Oyarte Galvez et al. (2025): experimental setup, tip speed 240 µm/h
