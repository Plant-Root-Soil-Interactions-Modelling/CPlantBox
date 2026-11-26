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

# ---- test: nominal point over 180 h (LONG_HOURS), roots turned off at t* and roots growing on (10 seeds each) ----
# are the wave fits defined often enough, when do new hyphae stop coming without root growth, from when on do the outputs differ?
sbatch sa_run.sbatch design.py long          # wait until the job has finished
SIMTIME=24:00:00 bash submit_sa.sh results/design_long.csv   # roots on is probably slower: possibly more time
# result: slurmOut/sa_collect-<job>.out (analyze_long.py), windows 60, 90 and 120 h, results/long_roots_compare.csv,
# long_influx.csv and plots/long_influx.png (supply of hyphae per hour). Long runs from before the supply diagnostics
# (raw.npz without diag_*) are reported; delete runs/p0_long_*/ and resubmit to include them

# ---- Phase 1: Morris (r × (k+1) × R runs) ----
sbatch sa_run.sbatch design.py morris        # wait until the job has finished
bash submit_sa.sh results/design_morris.csv
# analyze_morris.py runs automatically at the end, plots in sa_setup/results/plots/

# ---- inoculum Morris: root phase, root architecture, placement (r × 6 × R runs, hyphal parameters nominal) ----
sbatch sa_run.sbatch design.py inoculum      # wait until the job has finished
bash submit_sa.sh results/design_inoculum.csv
# analyze_morris.py inoculum runs automatically: results/inoculum_indices.csv, plots/inoculum_*.png and, once the
# hyphal Morris exists, inoculum_vs_hyphae.csv (max μ* root/inoculum / max μ* hyphae per output)
```

**Inoculum Morris** (`config.INOCULUM_FACTORS`): the root phase (`ROOT_DAYS`), the root architecture (tomato as a
stand-in for the root organ culture of *D. carota*, no data on its geometry; factors on `r` and `ln` of all roots) and
the placement of the roots (`seed_pos_x`) have no source. Plus the hyphal emergence density (Schnepf et al.) for comparison.
Question: does the normalisation to t* remove their effect? They may change `wait_days` and the supply of hyphae
(`influx_tstar`, `influx_window`, `influx_decay` from the crossings per hour), but should hardly change the outputs after t*.
Runs with roots turned off at t*, as SA and calibration.

**Ring resolution** (`analyze_rings.py`, no re-simulation): merges every 3 neighbouring rings (15 -> 5, again equal
area) and compares all ring-dependent outputs, per run (SA) and pooled (calibration targets):
`sbatch sa_run.sbatch analyze_rings.py` (Phase 0 noise runs) or `sbatch sa_run.sbatch analyze_rings.py results/design_long.csv long 120`
-> `results/rings_compare_<design>.csv`.

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

## Parameters (one list for SA and calibration)

`config.PARAMETERS` defines every parameter once: range (= Morris range = calibration search range), prior, source and
`calibrate` (free in the calibration or fixed at the nominal value). The SA varies all of them plus `config.SA_ONLY`.

| Group | Parameters | Range |
|---|---|---|
| Runner (subType 1) | `R_v`, `R_hlt`, `R_lmax`, `R_ln` (lb = ln), `R_theta`, `R_tropismN`, `R_succ_BAS`, `R_a` | v 0.432–0.72 cm/d, hlt 4–10 d, lmax 1.5–5 cm, ln 0.05–0.6 cm (log), θ 50–80°, tropismN 1–8, BAS share of the laterals 0.5–0.95, radius 10–100 µm (log) |
| BAS (subType 2) | `B_v`, `B_b`, `B_hlt`, `B_theta`, `B_tropismN`, `B_a`, `B_lmax` | v as runner, b 0.3–2 /d (log), hlt 3–12 d, θ 50–80°, tropismN 0.5–4, radius 0.5–50 µm (log), lmax 0.01–0.1 cm (log) |
| shared | `ana`, `distTH`, `ana_L_scale` | 0.5–1, 0.005–0.03 cm (log), anaMinDist factor 0.5–2 on the derived value (log) |
| SA only | `v_cv`, `cv_scale`, `R_b_ctrl` | CV of v 0.15–0.40, all other dev × 0.5–2 (log), control (b of the runners, no effect with b_prob = 0) |

21 factors: 20 trajectories × 22 points × 5 replicates = 2200 runs. Fixed: la, lb of the BAS, tropismT, tropismS, dx, b_prob.

**First Morris round:** its 13 factors are kept as `config.FACTORS_V1`. Designs without the column `factor_set` (first
round) are run and evaluated with them, so its runs stay usable (e.g. `--eval-only` for the pooled targets). New designs
store `factor_set = v2` and have new run_ids (`p0v3_*`, `m2_*`, `inoc3_*`): the nominal model changed (anaMinDist
derived; 120 h after t*).

**Order (B14):** Phase 0 (noise, dx, dt, layer, anaMinDist) and the long test -> fix the numerics; Morris (and the
inoculum Morris) -> `analyze_morris.py` writes `morris_calibration_suggestion.csv` (parameters with μ* above noise for
a calibration target). **Only a suggestion:** set `calibrate=True/False` in `config.PARAMETERS` by hand, then calibrate.
The fixed parameters are listed in `calibration/report.txt`.

## Rings

As in Oyarte Galvez et al. 2025 the radial densities are computed against the distance r from the **barrier crossing
point of each run** (median of the hyphal crossings of the barrier until t*, `AMFSimulationSetup.crossingPoint`), in
`N_RINGS` = 15 rings of equal area up to `RING_RADIUS` = 3.2 cm (fits the ring positions r ≈ 11–31 mm of Fig. 2c; derived,
not stated). Areas are computed numerically (`ringsAround`): off centre the dish edge can cut the outer rings.
`RING_CENTER = "dish"` restores the rings around the centre of the dish. `analyze_phase0.py` reports where the crossing
point lies.

**Targets from a single colony:** Fig. 2 shows one colony, its bands are spatial bootstraps within a ring (10,000
rectangles), not the spread between plates. These targets (`replicates=1`: K1, K2, `front_velocity`, `front_tip_speed`)
get `between_cv`, the CV between the seeds at the nominal point that `analyze_phase0.py` prints, entered by hand after
Phase 0: effective tolerance √(tol² + between_cv²).

## Outputs

| Kind | Outputs | Time |
|---|---|---|
| quantities, per run | `total_length` [cm], `active_tips`, `bas_share`, `anastomoses` (since t*), `front_position` (tip-weighted radius, mm), `spread_r90`, `max_tip_density` | 10, 20, 30, 40, 60, 90, 120 h after t* (`SNAPSHOT_HOURS`) |
| wave, per run | `front_speed` [mm/d], `encounter_C` | 0–120 h, stored at 120 h |
| calibration targets, pooled per design point | `alpha` [µm/h], `beta` [h⁻¹], `K1` [µm mm⁻²], `K2` [mm⁻²], `front_velocity` [mm/d], `wave_width` [h] | 0–120 h, stored at 120 h |
| image graph, per run | `root_efficiency_rel`, `bas_share_image`, `L_RH_image`, `L_BAS_image` [mm], `bas_share_subtype_image`, `ana_per_length_rh`, `ana_per_length_total` | 20, 40, 60 h after t* (`IMAGE_HOURS`) |

The image graph quantities are computed in `model_hook` from the final plant with the segments created until 20, 40 and
60 h after t* (planar graph as the image analysis, `AMFNetwork.getPlanarNetwork(t_max=...)`). The calibration uses their
time courses (and `anastomoses`) as targets at these times, read off the figures of Oyarte Galvez et al. 2025.

**One definition for SA and calibration** (`evaluate.WAVE_QUANTITIES`, `evaluate.wave_quantities`): tips are the
growing runner tips (as in the experiment), `front_velocity` from the arrival t_n of the filament density (sigmoid fit).
**Nomenclature:** `alpha` = anastomosis (a = α ρ n, target ≈ 23 µm/h), `beta` = branching (new runner tips per growing
runner tip and hour, target ≈ 4 %/h). Oyarte Galvez et al. (2025) use the letters the other way round.

The wave fits come from `AMFDensities` with the quality criteria of `petridishstats.py` (R², saturation, at least
3 rings). Per run they fail too often for Morris, so the calibration targets are **pooled per design point** over its
replicates, as the calibration pools over its seeds (`evaluate_pooled.py` -> `results/pooled_<design>.csv`, runs
automatically in `collect.sbatch`; noise level from a bootstrap over the replicates). `bas_share` is a target too
(snapshot at 120 h, per run).

**Simulated time (`HOURS_HYPHAE` = 120 h after t*):** the Phase 0 long test showed that K1 and K2 only reach the
experimental values at about 120 h (nominal K1 763 -> 1422 µm mm⁻², K2 0.071 -> 0.151 mm⁻² between 60 and 120 h), so the
wave fits run over 0–120 h. The image graph targets stay at 20/40/60 h (time course of the figures). The long test runs
`LONG_HOURS` = 180 h (image graph additionally at 150/180 h). Run ids changed to `p0v3_*` and `inoc3_*`; Phase 0 has to
be rerun with 120 h, the runtime per run is about 1.6 h (nominal median of the long test). Per run the wave quantities are still in `outputs.csv` as diagnostics, with the number of rings
behind the fits (`n_wave_rings`, `n_K1_rings`, `n_tip_rings`) and `wait_days`.

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
  checked once per step at the tip, `ana` and the BAS tip splitting act per step: check `alpha`, `anastomoses`,
  `encounter_C` and `bas_share`. Small dt costs runtime roughly ∝ steps: increase `--time` in `sim.sbatch`.
- `phase0_layer.csv`: thickness of the hyphal layer on the fungal side (0.1 / 0.05 / 0.02 cm,
  `config.PHASE0_LAYER_VALUES`) vs. the thinnest, closest to the quasi-2D cellophane. Check `alpha`, `anastomoses`,
  `encounter_C`; aborts in thin layers show up in `phase0_status.csv`.
- `phase0_anaL.csv`: minimal path length along the mycelium for anastomosis (`anaMinDist`). Nominal now derived by the
  model (`-1` in the template: distTH / sin θ + 2 dx ≈ 0.03 cm instead of the former default 0.1 cm), tested as a factor
  0.5 / 2 / 3.17 (≈ 0.1 cm) on the derived value. For the next Morris round a factor with param `ana_L_scale` is supported.
- `frozen_share` (all designs, from the diagnostics): share of the tips on the fungal side that are frozen because an
  ancestor reached its lifetime `hlt` (`Hyphae::simulate` simulates the children of a hypha only while it is alive).
  Grows with the waiting time until t*.
- `phase0_crn.csv`: common random numbers. `PHASE0_CRN_PARAM` moved by the calibration step (0.15 normalised), once with
  the seeds of the noise runs, once with independent seeds: sd of the difference same / independent seeds. Clearly < 1:
  shared seeds remove noise in the calibration differences and the Morris steps; ~1: the runs decorrelate, more seeds needed.
- `phase0_bas_share.csv`: BAS share by sub type vs. as the experiment classifies the image graph (C15).
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

Parameters: those of `config.PARAMETERS` with `calibrate=True` (all 18 by default, set by hand after the Morris
screening; the others stay at the nominal value). Separate for runner/BAS (v, hlt, theta, ln, lmax, tropismN, radius a,
BAS share of the runner successors, b and lmax of the BAS: a BAS segment splits only while it grows, so `b · lmax / v`
sets the size of the BAS), shared: ana, distTH, `ana_L_scale`. lb of the runners is coupled to ln, la stays fixed.
Targets and tolerances are read off the figures of Oyarte Galvez et al. 2025 (`Target.source`; approximate). t* is, as
in the paper, the time at which the runner length behind the barrier reaches 10 cm (`crossing_subtype = 1` in
`makeSimulation`). Additional targets: `L_RH_image` ≈ 800 mm and `L_BAS_image` ≈ 450 mm at 60 h after t* (lengths as
the image analysis classifies them). The anastomosis density ρ_A is computed in two variants (`ana_per_length_rh`,
`ana_per_length_total`) until its definition in the paper is checked.
Current targets: `alpha` (anastomosis) 23 µm/h, K1 1400 µm mm⁻², BAS share 0.38 (`bas_share_image`: classified as in
the experiment from the edges of the planar image graph, edge < 400 µm, or < 1000 µm ending in a tip, or width < 7 µm
and width × length < 9000 µm², `AMFNetwork.basShareImage`), K2 0.15 mm⁻² (runner tips), `beta`
(branching) 4 %/h (new runner tips per growing runner tip and hour, `branch_new_runner` in `raw.npz`),
wave speed 280 µm/h (all Oyarte Galvez et al. 2025) and relative root efficiency
≈ 0.5: (E_r − E_r,MST)/(E_r,DT − E_r,MST) between tips/BAS attachment sites and root nodes (≤ 5 mm
from the barrier), at the end of the simulation (`AMFNetwork.relativeRootEfficiency`, can be switched off via
`config.ROOT_EFFICIENCY`). As in the experiment, a planar graph seen from above is evaluated
(`AMFNetwork.getPlanarNetwork`): all hyphal segments on the fungal side projected onto the x-y plane, every
crossing of two segments becomes a node (degree 4 in the image skeleton).
Because 7 targets do not determine 18 parameters uniquely, every parameter has a **prior**, an
additional residual z = (x − centre)/σ (log for log parameters). Literature priors (`Factor.prior`, range
= ±2σ) exist for R_v, B_v, R_hlt, R_ln, R_theta, B_theta, B_hlt, B_a (BAS diameter 2–10 µm, Bago 1998) and R_a
(radius 50 µm, Schnepf 2017, prior 25–100 µm; lower bound of the search range 10 µm, Bago 1998), and an assumed one for
`ana_L_scale` (0.5–2 around the derived anaMinDist). All others have a weak
prior around the nominal value (σ = `WEAK_PRIOR_SIGMA` × search range). Parameters that no target acts on
thus stay at the nominal value. The report flags literature priors from which the solution deviates by more than 2σ.
All evaluations use the same seeds (CRN), the outputs are pooled over the seeds, the
standard error comes from a bootstrap. Parameters, search ranges and priors are in `config.PARAMETERS`, targets and
tolerances at the top of `calibrate.py` (`TARGETS`, with `source`: most tolerances are still assumptions, they should be
the experimental uncertainty of the target).

`report.txt` also lists per parameter the **data share** (information from the targets / (targets + prior): near 0 means
the "calibrated" value is its prior) and a **±1σ range** at the best point (linearised, from the last Jacobian; σ is the
change that moves the targets by about one tolerance). By default the search stops at the first point with all targets
within tolerance; `--full` goes on to the minimum of the cost.

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
