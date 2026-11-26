"""
Calibration of the hyphal parameters to target values as a root-finding problem, Petri dish setup as in the SA.

    python calibrate.py --check            # print parameters, targets, start point, cost
    python calibrate.py                    # calibration (under SLURM: sbatch calibrate.sbatch)
    python calibrate.py --seeds 8 --max-nfev 10

We look for the parameter vector x with F(x) = 0, where

    F_i(x) = log(y_i(x) / y*_i) / log(1 + tol_i)

y_i: output i, pooled over N_SEEDS seeds, y*_i: target value, tol_i: accepted relative deviation.
|F_i| <= 1 means "within tolerance". The parameters are searched in the normalised space u in [0,1]
(log-uniform as in design.to_physical), so that steps are comparable between parameters.

Solver: scipy.optimize.least_squares (Trust Region Reflective, box bounds [0,1]). With at least as many
targets as parameters, the minimum of |F|^2 is the root; scipy.optimize.root does not support bounds.
With fewer targets than parameters the solution is not unique: the script then stops
(--allow-underdetermined forces the run). The Jacobian comes from forward differences with a large
step DIFF_STEP (the signal must be clearly above the seed noise); all difference points
are distributed over the cores together. The search stops early as soon as all |F_i| <= 1.

Noise: all evaluations use the same seeds (common random numbers), so differences
between parameter points are less noisy. The standard error of each output comes from a bootstrap
over the seeds: if |y - y*| is below ~2 SE, the deviation is within the noise.

Output in calibration/: evals.csv (every evaluation), jacobian.csv (sensitivities dF_i/du_j, columns
close to 0 = parameter not determined by the targets), best_params.xml (best point). Restartable: finished
runs (runs/<point>/seed_<s>/meta.json) are not simulated again, the solver replays deterministically.
"""
import argparse
import hashlib
import json
import math
import os
import signal
import time
import traceback
import warnings
from collections import Counter
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

import config
import evaluate as ev
import model_hook
from xml_tools import factor_nominal, write_xml

amfd = ev.amfd


# --------------------------------------------------------------------------
# Parameters: the shared list config.PARAMETERS (same ranges as the SA, B13). Free are those with calibrate=True,
# chosen BY HAND after the Morris screening (analyze_morris.py writes a suggestion); the others stay at the nominal value.
# --------------------------------------------------------------------------
PARAMS = [p for p in config.PARAMETERS if p.calibrate]
FIXED_PARAMS = [p for p in config.PARAMETERS if not p.calibrate]
# Weak prior for parameters without literature: centre = nominal value, σ = WEAK_PRIOR_SIGMA × search range
# (in log for log parameters). Keeps them at the nominal value as long as the targets do not need them; None = no prior.
WEAK_PRIOR_SIGMA = 0.5


# --------------------------------------------------------------------------
# Targets
# --------------------------------------------------------------------------
@dataclass
class Target:
    output: str                 # see OUTPUTS
    value: float
    tol_rel: float              # accepted relative deviation
    hour: float | None = None   # snapshot outputs only: time [h after t*]
    log: bool = True            # F = log(y/y*)/log(1+tol); False: F = (y - y*)/(tol y*), for outputs that can be <= 0
    # where the value and the tolerance come from. The tolerance should be the experimental uncertainty of the target
    # (e.g. spread between plates or 2 SEM, relative): then |F| <= 1 means "within the experimental uncertainty" and the
    # parameter uncertainty in report.txt (C18) refers to it. "assumption" = set by hand (C16)
    source: str = "assumption"
    # tolerance widened by this factor in the residual (F / widen): for the time courses √(number of time points), since
    # the time points of one replicate are correlated and together should weigh about as much as one target
    widen: float = 1.0
    # "lower": the value is a lower bound (F = 0 above it), e.g. the persistence length of the paper
    bound: str | None = None
    # number of replicates (plates) behind the value: 1 = a single colony (Fig. 2 of the paper: its bands are spatial bootstraps
    # within one colony, not the spread between plates); None = unknown
    replicates: int | None = None
    # spread between plates (relative), for targets from a single colony, which the experiment does not give: estimated from
    # the model, CV of the output between the seeds at the nominal point (phase0_noise.csv), set by hand after Phase 0.
    # Effective tolerance sqrt(tol_rel² + between_cv²)
    between_cv: float = 0.0

    @property
    def tol_eff(self):
        return float(np.hypot(self.tol_rel, self.between_cv))

    def residual(self, y):
        if self.bound == "lower" and y == np.inf:
            return 0.
        if not np.isfinite(y) or (self.log and y <= 0):
            return np.nan
        if self.bound == "lower" and y >= self.value:
            return 0.
        tol = self.tol_eff
        F = np.log(y / self.value) / np.log(1 + tol) if self.log else (y - self.value) / (tol * self.value)
        return F / self.widen

    @property
    def tol_text(self):
        return f"{self.tol_rel:.0%}" + (f" ⊕ CV {self.between_cv:.0%}" if self.between_cv else "") + \
            (f" ×{self.widen:.2f}" if self.widen != 1 else "") + (" (lower bound)" if self.bound == "lower" else "") + \
            (" [1 colony]" if self.replicates == 1 else "")

    @property
    def label(self):
        return self.output if self.hour is None else f"{self.output}@{self.hour:g}h"


# Available outputs. Wave quantities over the window 0..HOURS, pooled over all seeds, same definitions as in the SA
# (evaluate.WAVE_QUANTITIES). Nomenclature: alpha = anastomosis, beta = branching (Oyarte Galvez et al. 2025: the other way round)
WAVE_OUTPUTS = ev.WAVE_QUANTITIES
# Snapshot outputs (mean over seeds) at time Target.hour: config.SNAPSHOT_OUTPUTS
# total_length [cm], active_tips, bas_share, anastomoses (since t*), front_position [mm], spread_r90 [mm],
# max_tip_density [mm^-2]; at config.IMAGE_HOURS also the image graph quantities (evaluate.IMAGE_QUANTITIES)

# A unique solution needs at least as many targets as parameters (here 18), and every
# parameter must clearly affect at least one target. Add further measurements as Target(...).
TS = math.sqrt(len(config.IMAGE_HOURS))   # widening of the tolerances of the time courses (correlated time points)
TARGETS = [
    # a(t) and b(t) against t − t_n (Fig. 3b, scaled onto n and ρn): the bands of a are about ±40–50 % per time point,
    # of b about ±30–40 %; the constant is fitted over all points, so somewhat narrower. Fig. 3 shows ONE replicate:
    # replicates=1, between_cv (spread between plates, from the seeds in Phase 0) is added after Phase 0
    Target("alpha", 23.0, 0.30, replicates=1, between_cv=0.0,
           source="Fig. 3b: anastomosis a = α ρ n, ≈ 23 µm/h (their beta); band of a(t) ≈ ±40–50 % per time point (one replicate)"),
    Target("beta", 0.04, 0.30, replicates=1, between_cv=0.0,
           source="Fig. 3b: branching b = β n, ≈ 4 %/h (their alpha); band of b(t) ≈ ±30–40 % per time point (one replicate)"),
    # values read off the figures of Oyarte Galvez et al. 2025 (approximate, check against the data): Fig. 2 and the figure
    # with L_RH, L_BAS, ρ_BAS, anastomoses and Ê_g, Ê_r, Ĉ over t − t* (12 replicates, mean ± s.d.)
    # Fig. 2 shows ONE colony; its bands are spatial bootstraps within a ring (10,000 rectangles) or fit intervals, not the
    # spread between plates. replicates=1; between_cv (spread between plates, from the model's seeds at the nominal point,
    # phase0_noise.csv) is set by hand after Phase 0 — until then 0, i.e. these targets are too strict
    Target("K1", 1400.0, 0.12, replicates=1, between_cv=0.0,
           source="Fig. 2c: saturation density ≈ 1400 µm mm⁻², points/band ≈ ±150–200 between rings (one colony)"),
    Target("K2", 0.15, 0.33, replicates=1, between_cv=0.0, source="Fig. 2b: tip peak 0.1–0.2 mm⁻² (one colony)"),
    Target("front_velocity", 280e-3 * 24, 0.10, replicates=1, between_cv=0.0,
           source="Fig. 2d inset: wave front r vs. t_n ≈ 280 µm/h = 6.72 mm/d (one fit, one colony); tolerance: assumption"),
    # growth speed of the tips at the front (outermost 10 % of the growing runner tips, model_hook.tip_speed_diagnostics),
    # mean from 10 h after t*; together with the prior on v (mean ≈ 240) it determines v_cv
    Target("front_tip_speed", 300.0, 0.10, replicates=1, between_cv=0.0,
           source="Fig. 2d: average growth speed at the front ≈ 300 µm/h, constant; points per time step ≈ 250–330 (one colony)"),
    # Extended Data Fig. 11. Angles: the mean of 881 points is precise, the tolerance covers the unknown way of measuring
    # (projection, length over which the directions are taken); the priors on theta are weak, so the data count once
    Target("branch_angle_rh", 73.0, 0.10, source="Ext. Data Fig. 11b: new-born RH branch at 73° from the initial growth direction"),
    Target("branch_angle_bas", 80.0, 0.10, source="Ext. Data Fig. 11b: new-born BAS at 80° (not significantly different from RH)"),
    Target("persistence_length", 10.0, 0.5, bound="lower",
           source="Ext. Data Fig. 11a: median persistence length ≈ 10 cm, a lower bound (decay barely visible up to 1 cm)"),
    # ---- time courses at 20, 40, 60 h after t* (config.IMAGE_HOURS), mean of 12 replicates, tolerance from the s.d. band ----
    # the time points of one replicate are correlated: each tolerance is widened by √3 (TS), so that the three points of a
    # quantity together weigh about as much as one target
    # runner length as the image analysis classifies it (BAS criteria of the edges, projected length) [mm], log axis
    Target("L_RH_image", 280.0, 0.30, hour=20, widen=TS, replicates=12, source="L_RH ≈ 280 mm at 20 h, band ≈ ±30 %"),
    Target("L_RH_image", 500.0, 0.30, hour=40, widen=TS, replicates=12, source="L_RH ≈ 500 mm at 40 h, band ≈ ±30 %"),
    Target("L_RH_image", 850.0, 0.40, hour=60, widen=TS, replicates=12, source="L_RH ≈ 850 mm at 60 h, band ≈ ±40 %"),
    # BAS share as classified in the experiment (edge < 400 µm, or < 1000 µm ending in a tip, or width < 7 µm and
    # width × length < 9000 µm²); L_BAS (≈ 100, 230, 370 mm) follows from L_RH and ρ_BAS and is only an output
    Target("bas_share_image", 0.37, 0.14, hour=20, widen=TS, replicates=12, source="ρ_BAS ≈ 0.37 at 20 h, s.d. ≈ 0.05"),
    Target("bas_share_image", 0.385, 0.13, hour=40, widen=TS, replicates=12, source="ρ_BAS ≈ 0.385 at 40 h, s.d. ≈ 0.05"),
    Target("bas_share_image", 0.40, 0.13, hour=60, widen=TS, replicates=12, source="ρ_BAS ≈ 0.40 at 60 h, s.d. ≈ 0.05"),
    # cumulative anastomoses on the fungal side; the model counts from t*, the experiment has ≈ 21 at t*, subtracted; log axis
    Target("anastomoses", 75.0, 0.50, hour=20, widen=TS, replicates=12, source="anastomoses ≈ 95 at 20 h (≈ 21 at t*), band ≈ ±35–50 %"),
    Target("anastomoses", 240.0, 0.60, hour=40, widen=TS, replicates=12, source="anastomoses ≈ 260 at 40 h (≈ 21 at t*), band ≈ ±60 %"),
    Target("anastomoses", 445.0, 0.60, hour=60, widen=TS, replicates=12, source="anastomoses ≈ 465 at 60 h (≈ 21 at t*), band ≈ ±55–60 %"),
    # relative root efficiency (can be < 0: linear residual, tol_rel × value = absolute tolerance)
    Target("root_efficiency_rel", 0.69, 0.15, hour=20, widen=TS, replicates=12, log=False, source="Ê_r ≈ 0.69 at 20 h, s.d. ≈ 0.1"),
    Target("root_efficiency_rel", 0.665, 0.12, hour=40, widen=TS, replicates=12, log=False, source="Ê_r ≈ 0.665 at 40 h, s.d. ≈ 0.08"),
    Target("root_efficiency_rel", 0.63, 0.13, hour=60, widen=TS, replicates=12, log=False, source="Ê_r ≈ 0.63 at 60 h, s.d. ≈ 0.08"),
    # not targets: L_BAS (redundant), anastomosis density ρ_A ≈ 0.53, 0.65, 0.69 (the paper does not define it; besides
    # redundant with anastomoses and L_RH), global efficiency Ê_g ≈ 0.33, 0.44, 0.53
]

HOURS = config.HOURS_HYPHAE   # simulation time after t* [h], window of the wave quantities
N_SEEDS = 16                  # seeds per parameter point
SEED0 = config.BASE_SEED + 300_000
DIFF_STEP = 0.15              # difference step in the normalised space
MAX_NFEV = 30                 # max. solver iterations (least_squares: max_nfev, without differences)
MIN_OK = 0.5                  # share of seeds that must reach t*, otherwise the point counts as invalid
NAN_PENALTY = 10.             # F_i for undefined outputs (fit failed), so that the solver keeps going
# True: stop at the first point with all |F| <= 1 (cheap). False (--full): go on to the minimum of |F|² + |z|², so that
# the result does not depend on where the search entered the tolerance region and the ±1σ in report.txt refer to the minimum
STOP_AT_TOLERANCE = True
N_BOOT = 200
CAL_DIR = config.ROOT / "calibration"
# Version of the raw data (raw.npz) in the key of the runs: increase when model_hook stores new quantities that the
# targets need; old runs are then not reused. 2: branch_new_runner, 3: bas_image_share (BAS share of the image graph),
# 4: t* from the runner length only (as L_RH = 10^2 mm in Oyarte Galvez et al. 2025), 5: image graph at IMAGE_HOURS, tip speeds,
# 6: rings around the barrier crossing point up to RING_RADIUS
RAW_VERSION = 6


# --------------------------------------------------------------------------
# Parameter space
# --------------------------------------------------------------------------
def spec_nominal(p):
    return factor_nominal(p)


def prior_spec(p):
    """(centre, σ, kind) on the scale of the parameter (log or linear), or None."""
    g = np.log if p.log else (lambda v: v)
    if p.prior is not None:
        lo, hi = p.prior
        return (g(lo) + g(hi)) / 2, (g(hi) - g(lo)) / 4, "literature"
    if WEAK_PRIOR_SIGMA:
        return g(spec_nominal(p)), WEAK_PRIOR_SIGMA * (g(p.high) - g(p.low)), "weak"
    return None


def prior_residuals(u):
    """Additional residuals (g(x) - centre) / σ per parameter with a prior: |z| = 1 means one standard deviation."""
    z = []
    for p, xj in zip(PARAMS, to_phys(u)):
        s = prior_spec(p)
        if s:
            z.append(((np.log(xj) if p.log else xj) - s[0]) / s[1])
    return np.array(z)


def prior_jacobian():
    """dz/du is constant: g(x) is linear in u."""
    rows = []
    for j, p in enumerate(PARAMS):
        s = prior_spec(p)
        if s:
            row = np.zeros(len(PARAMS))
            row[j] = (np.log(p.high / p.low) if p.log else p.high - p.low) / s[1]
            rows.append(row)
    return np.array(rows).reshape(-1, len(PARAMS))


def to_phys(u):
    return [float(np.exp(np.log(p.low) + uj * np.log(p.high / p.low)) if p.log else p.low + uj * (p.high - p.low))
            for p, uj in zip(PARAMS, u)]


def to_unit(x):
    return np.array([np.log(xj / p.low) / np.log(p.high / p.low) if p.log else (xj - p.low) / (p.high - p.low)
                     for p, xj in zip(PARAMS, x)])


def xml_updates(x):
    """Physical parameter values -> {parameter name: value} for xml_tools.write_xml (fixed parameters: template value)."""
    return {p.name: float(xj) for p, xj in zip(PARAMS, x)}


def write_params(x, path):
    return write_xml(xml_updates(x), path, factors=config.PARAMETERS)


# --------------------------------------------------------------------------
# Simulation and outputs
# --------------------------------------------------------------------------
def last_line(tb):
    lines = [l for l in str(tb).strip().splitlines() if l.strip()]
    return lines[-1].strip()[:300] if lines else "?"


def run_seed(args):
    """One run (in the worker process). Returns (path to raw.npz or None, error message or None)."""
    out = lambda meta: ((str(rdir / meta["snapshots"][0]), None) if meta["sim_status"] == "ok"
                        else (None, last_line(meta.get("error", "?"))))
    updates, seed, rdir, hours = args
    meta_path = rdir / "meta.json"
    if meta_path.exists():
        return out(json.loads(meta_path.read_text()))
    xml_path = write_xml(updates, rdir / "params.xml", factors=config.PARAMETERS)
    meta = {"seed": seed, "hours": hours, "values": updates}
    t0 = time.time()
    try:
        snaps = model_hook.simulate(xml_path, seed, rdir, hours=hours)
        meta.update(sim_status="ok", snapshots=[p.name for p in snaps])
    except Exception:
        meta.update(sim_status="failed", error=traceback.format_exc())
    meta["sim_seconds"] = round(time.time() - t0, 1)
    meta_path.write_text(json.dumps(meta, indent=1))
    return out(meta)


def prepare(sim):
    """Expensive evaluations once per run (fits, snapshots), the bootstrap reuses them.
    The wave fits use the density of the runner tips: in the experiment (Oyarte Galvez et al. 2025, Methods) only
    growing tips count (moved > 40 pixels), BAS tips hardly; Fig. 2/3 show the density of the growing RH tips."""
    hours = sorted({t.hour for t in TARGETS if t.hour is not None})
    return dict(ev.prepare(sim, HOURS), snap={h: ev.snapshot_outputs(sim, h) for h in hours})


def pooled_outputs(preps):
    wave = ev.wave_quantities(preps)   # same definitions as the SA (evaluate.WAVE_QUANTITIES)
    out = {}
    for t in TARGETS:
        if t.hour is None:
            out[t.label] = float(wave[t.output])
        else:
            out[t.label] = float(np.nanmean([p["snap"][t.hour].get(t.output, np.nan) for p in preps]))
    return out


def _outside_at(sim, hour):
    """Share of the growing runner tips beyond RING_RADIUS at `hour` after t* (model_hook.tip_speed_diagnostics), nan without."""
    if "diag_frac_tips_outside_rings" not in sim:
        return np.nan
    h = np.asarray(sim["diag_speed_hours"])
    return float(np.asarray(sim["diag_frac_tips_outside_rings"], dtype=float)[int(np.argmin(np.abs(h - hour)))])


class Converged(Exception):
    pass


class Residual:
    """F(u) and Jacobian for least_squares. All runs of one call (one point or all
    difference points) go to the pool together."""

    def __init__(self, pool, n_seeds, hours):
        self.pool, self.seeds, self.hours = pool, [SEED0 + i for i in range(n_seeds)], hours
        self.records, self.cache = [], {}
        self.best = None    # (|F|, u, F)
        self.errors, self.eff_errors = Counter(), Counter()   # error messages of the runs and of the root efficiency
        self.n_runs, self.J, self.status, self.t_start = 0, None, "running", time.time()

    def key(self, u):
        return hashlib.md5(json.dumps([round(v, 10) for v in to_phys(u)] + [self.hours, RAW_VERSION]).encode()).hexdigest()[:10]

    def evaluate(self, us, kind):
        """F for several points; already computed ones from the cache."""
        todo = [u for u in us if self.key(u) not in self.cache]
        jobs, owner = [], []
        for i, u in enumerate(todo):
            upd = xml_updates(to_phys(u))
            pdir = CAL_DIR / "runs" / self.key(u)
            jobs += [(upd, s, pdir / f"seed_{s}", self.hours) for s in self.seeds]
            owner += [i] * len(self.seeds)
        t0 = time.time()
        results = self.pool.map(run_seed, jobs, chunksize=1)
        minutes = round((time.time() - t0) / 60, 1)
        self.n_runs += len(results)
        self.errors.update(e for _, e in results if e)
        for i, u in enumerate(todo):
            mine = [p for (p, _), o in zip(results, owner) if o == i and p]
            self.cache[self.key(u)] = self.score(u, mine, kind, minutes)
        return [self.cache[self.key(u)] for u in us]

    def score(self, u, paths, kind, minutes):
        labels = [t.label for t in TARGETS]
        preps = [prepare(ev.truncate(model_hook.load_raw(Path(p)), self.hours)) for p in paths]
        self.eff_errors.update(last_line(p["sim"]["root_eff_error"]) for p in preps if "root_eff_error" in p["sim"])
        ok = len(preps) >= MIN_OK * len(self.seeds)
        y = pooled_outputs(preps) if ok else dict.fromkeys(labels, np.nan)
        rng = np.random.default_rng(int(self.key(u), 16) % 2 ** 32)
        boot = [pooled_outputs([preps[i] for i in rng.integers(0, len(preps), len(preps))])
                for _ in range(N_BOOT if ok else 0)]
        se = {n: float(np.nanstd([b[n] for b in boot])) if boot else np.nan for n in labels}
        F = np.array([t.residual(y[t.label]) for t in TARGETS])

        x = to_phys(u)
        z = prior_residuals(u)
        cost = float(np.sum(F ** 2) + np.sum(z ** 2))   # what least_squares minimises (× 2)
        rec = {"eval": len(self.records), "kind": kind, "point": self.key(u), "n_ok": len(preps), "minutes": minutes,
               **{p.name: xj for p, xj in zip(PARAMS, x)}, **{f"u_{p.name}": uj for p, uj in zip(PARAMS, u)},
               **y, **{f"se_{n}": se[n] for n in labels}, **{f"F_{n}": Fi for n, Fi in zip(labels, F)},
               "norm_F": float(np.linalg.norm(F)), "norm_prior": float(np.linalg.norm(z)), "cost": cost,
               # ring diagnostics for the warnings in report.txt (crossings until t*, share of tips beyond RING_RADIUS)
               "crossings_median": float(np.nanmedian([float(p["sim"].get("n_crossings_tstar", np.nan)) for p in preps]))
               if preps else np.nan,
               "outside_rings_median": float(np.nanmedian([_outside_at(p["sim"], self.hours) for p in preps])) if preps else np.nan}
        self.records.append(rec)
        pd.DataFrame(self.records).to_csv(CAL_DIR / "evals.csv", index=False)
        report(rec)
        if np.all(np.isfinite(F)) and (self.best is None or cost < self.best[0]):
            self.best = (cost, np.array(u), F)
        write_report(self)
        return F

    def targets(self, u):
        """Residuals of the targets (undefined -> NAN_PENALTY); convergence when all targets are within tolerance."""
        F = self.evaluate([u], "point")[0]
        if not np.all(np.isfinite(F)):
            if len(self.records) == 1:
                raise SystemExit("Not all outputs are defined at the start point (fits failed or too many runs "
                                 "without t*). Try more seeds (--seeds) or a longer simulation (--hours).")
            return np.where(np.isfinite(F), F, NAN_PENALTY)
        if STOP_AT_TOLERANCE and np.all(np.abs(F) <= 1):
            self.best = (0., np.array(u), F)   # the result is this point, even if another one had a lower cost
            raise Converged
        return F

    def __call__(self, u):
        """Targets and priors together: this is the residual vector for least_squares."""
        return np.concatenate([self.targets(u), prior_residuals(u)])

    def jac(self, u):
        """Targets: forward differences, all points in parallel, backward at the upper bound. Priors: analytic."""
        F0 = self.targets(u)
        h = np.where(u + DIFF_STEP <= 1, DIFF_STEP, -DIFF_STEP)
        pts = [u + h[j] * np.eye(len(u))[j] for j in range(len(u))]
        Fs = self.evaluate(pts, "diff")
        J = np.column_stack([(np.where(np.isfinite(Fj), Fj, F0) - F0) / h[j] for j, Fj in enumerate(Fs)])
        save_jacobian(J)
        self.J = J
        write_report(self)
        return np.vstack([J, prior_jacobian()])


def report(rec):
    pars = ", ".join(f"{p.name}={rec[p.name]:.3g}" for p in PARAMS)
    print(f"\n[{rec['eval']} {rec['kind']}] {pars}  ({rec['n_ok']}/{N_SEEDS} seeds ok)")
    for t in TARGETS:
        n = t.label
        flag = "ok" if abs(rec[f"F_{n}"]) <= 1 else ("within noise" if abs(rec[n] - t.value) <= 2 * rec[f"se_{n}"] else "")
        print(f"   {n:>18}: {rec[n]:10.4g} ± {rec[f'se_{n}']:.2g}  (target {t.value:g} ± {t.tol_text})  "
              f"F = {rec[f'F_{n}']:+.2f}  {flag}")
    print(f"   |F| = {rec['norm_F']:.3f}, |Prior| = {rec['norm_prior']:.3f}", flush=True)


HINTS = {
    "K1": "too few saturated rings with good fits -> longer simulation (--hours 180) or more seeds",
    "K2": "too few rings with a complete tip peak -> longer simulation (--hours 180) or more seeds",
    "wave_width": "as K2", "front_velocity": "as K2",
    "root_efficiency_rel": "root efficiency not computed -> see root efficiency errors",
    "alpha": "too few rings with >= ENCOUNTER_MIN_TIPS tips", "beta": "as alpha (or raw.npz without branch_new_runner)",
    "bas_share_image": "image graph not computed -> bas_image_error in raw.npz (needs config.ROOT_EFFICIENCY and networkx)",
    "L_RH_image": "as bas_share_image",
    "front_tip_speed": "no growing runner tips on the fungal side, or raw.npz without diag_front_tip_speed (speed_error)",
    "branch_angle_rh": "no branches on the fungal side or geometry_error in raw.npz", "branch_angle_bas": "as branch_angle_rh",
    "persistence_length": "no runners longer than PERSISTENCE_MIN_LENGTH on the fungal side or geometry_error",
}


def write_report(r, status=None):
    """calibration/report.txt: state, best point and error analysis, rewritten after every evaluation."""
    if status:
        r.status = status
    labels = [t.label for t in TARGETS]
    recs = pd.DataFrame(r.records)
    L = [f"Calibration, as of {time.strftime('%Y-%m-%d %H:%M')}, runtime {(time.time() - r.t_start) / 3600:.1f} h",
         f"Status: {r.status}"]
    if r.status == "running":
        L.append("  (If the job has ended and this still says 'running', it was killed: time limit or node failure. "
                 "Just resubmit, finished runs are reused.)")
    n_pts = int((recs["kind"] == "point").sum()) if len(recs) else 0
    n_prior = sum(1 for p in PARAMS if prior_spec(p))
    n_lit = sum(1 for p in PARAMS if p.prior is not None)
    L += [f"{len(PARAMS)} parameters, {len(TARGETS)} targets, {n_prior} priors ({n_lit} literature)"
          f"{' (UNDERDETERMINED)' if len(TARGETS) + n_prior < len(PARAMS) else ''}, {len(r.seeds)} seeds, {r.hours} h after t*",
          f"Evaluations {len(recs)} ({n_pts} iteration points, {len(recs) - n_pts} difference points), runs {r.n_runs}",
          "Fixed at the nominal value (config.PARAMETERS, calibrate=False): "
          + (", ".join(f"{p.name}={spec_nominal(p):.4g}" for p in FIXED_PARAMS) or "none"), ""]

    if r.best is not None:
        _, u, F = r.best
        b = recs.loc[(recs[[f"u_{p.name}" for p in PARAMS]] - u).abs().max(axis=1).idxmin()]
        L.append(f"Best point (evaluation {b['eval']}, |F| targets = {b['norm_F']:.2f}, |z| priors = {b['norm_prior']:.2f}):")
        z = iter(prior_residuals(u))
        for p, uj in zip(PARAMS, u):
            edge = "  <- AT BOUND" if uj <= 1e-6 or uj >= 1 - 1e-6 else ""
            s = prior_spec(p)
            pz = f"prior {s[2]:<10} z = {next(z):+5.2f}" if s else "no prior"
            L.append(f"  {p.name:>11} = {b[p.name]:<10.4g} [{p.low:.4g}, {p.high:.4g}], nominal {spec_nominal(p):<8.4g} {pz}{edge}")
        L.append("Targets at the best point:")
        for t, Fi in zip(TARGETS, F):
            y, se = b[t.label], b[f"se_{t.label}"]
            flag = "ok" if abs(Fi) <= 1 else ("within noise" if abs(y - t.value) <= 2 * se else "off")
            L.append(f"  {t.label:>20}: {y:10.4g} ± {se:<8.2g} target {t.value:g} ± {t.tol_text}   F ={Fi:+6.2f}  {flag}")
        L.append("History of iteration points (cost |F|² + |z|²): " + ", ".join(f"{v:.1f}" for v in recs.loc[recs.kind == "point", "cost"]))
        L.append(f"-> best_params.xml" if (CAL_DIR / "best_params.xml").exists() else "")
    else:
        L.append("No point with all outputs defined yet.")

    L += ["", "CHECKS"]
    n_fail = sum(r.errors.values())
    L.append(f"Failed runs: {n_fail} of {r.n_runs}" + (" -> most frequent errors:" if n_fail else ""))
    L += [f"  {c:4d} x {msg}" for msg, c in r.errors.most_common(5)]
    if len(recs):
        low = recs[recs.n_ok < MIN_OK * len(r.seeds)]
        if len(low):
            L.append(f"  {len(low)} evaluations with fewer than {MIN_OK:.0%} successful seeds (count as invalid)")
        for t in TARGETS:
            nan = int(recs[t.label].isna().sum())
            if nan:
                L.append(f"Undefined: {t.label} in {nan} of {len(recs)} evaluations ({HINTS.get(t.output, '')})")
        for t in TARGETS:
            se = recs[f"se_{t.label}"].median()
            if np.isfinite(se) and 2 * se > t.tol_rel * abs(t.value):
                L.append(f"Noise: {t.label} SE {se:.3g} > half the tolerance -> more seeds or a larger tolerance")
    if r.best is not None:
        zs = iter(prior_residuals(r.best[1]))
        for p in PARAMS:
            s = prior_spec(p)
            zj = next(zs) if s else 0.
            if s and s[2] == "literature" and abs(zj) > 2:
                L.append(f"Conflict with literature: {p.name} is {zj:+.1f} σ from its prior ({p.source}) -> the targets "
                         f"require a different value; check target, model structure or prior")
    if len(recs) and r.best is not None:   # rings at the best point (medians over its seeds)
        from analyze_phase0 import ring_warnings
        b = recs.loc[(recs[[f"u_{p.name}" for p in PARAMS]] - r.best[1]).abs().max(axis=1).idxmin()]
        L += ring_warnings([b.get("crossings_median", np.nan)], [b.get("outside_rings_median", np.nan)], " (best point)")
    if r.eff_errors:
        L.append(f"Root efficiency errors in {sum(r.eff_errors.values())} runs:")
        L += [f"  {c:4d} x {msg}" for msg, c in r.eff_errors.most_common(3)]
    if r.J is not None:
        J = r.J
        col = np.abs(J).max(axis=0)
        dead = [p.name for p, c in zip(PARAMS, col) if c < 0.3]
        if dead:
            L.append("No clear influence on the targets (max |dF/du| < 0.3), determined only by their prior: " + ", ".join(dead))
        for t, row in zip(TARGETS, J):
            top = np.argsort(-np.abs(row))[:3]
            weak = "  <- NO parameter has a clear effect" if np.abs(row).max() < 0.3 else ""
            L.append(f"  {t.label:>20} is driven by " + ", ".join(f"{PARAMS[k].name} ({row[k]:+.2f})" for k in top) + weak)
        sv = np.linalg.svd(J, compute_uv=False)
        L.append(f"Target Jacobian: rank {int((sv > 0.05 * sv[0]).sum())} of {min(J.shape)} "
                 f"(singular values {', '.join(f'{s:.2g}' for s in sv)}); the priors supply the rest")
        L += identifiability(J, r.best[1] if r.best is not None else None)
    if not any(l.startswith(("Failed", "Undef", "Noise", "Root", "No clear", "Conflict", "ATTENTION")) and not l.startswith("Failed runs: 0")
               for l in L[L.index("CHECKS"):]):
        L.append("No issues found.")
    (CAL_DIR / "report.txt").write_text("\n".join(L) + "\n", encoding="utf-8")


def identifiability(J, u_best):
    """C17/C18, from the last target Jacobian J (dF/du) and the priors (linearised, Gauss-Newton):
    - data share per parameter: information from the targets / (targets + prior) = (JᵀJ)_jj / ((JᵀJ)_jj + (PᵀP)_jj).
      Near 1: determined by the targets, near 0: only by its prior (the "calibrated" value is the prior).
    - ±1σ range at the best point from the covariance (JᵀJ + PᵀP)⁻¹ in the normalised space, converted to physical values.
      F is scaled with the tolerances: σ is the parameter change that moves the targets by about one tolerance, i.e. a
      parameter uncertainty in the sense of the experimental uncertainty only if the tolerances are that (C16).
    J comes from the last Jacobian point, which can differ from the best point (linearisation)."""
    P = prior_jacobian()
    JtJ, PtP = J.T @ J, P.T @ P
    data = np.diag(JtJ) / np.maximum(np.diag(JtJ) + np.diag(PtP), 1e-12)
    try:
        sd_u = np.sqrt(np.clip(np.diag(np.linalg.inv(JtJ + PtP)), 0, None))
    except np.linalg.LinAlgError:
        sd_u = np.full(len(PARAMS), np.nan)
    L = ["Identifiability (linearised, last Jacobian): data share = targets / (targets + prior), ±1σ at the best point:"]
    u = u_best if u_best is not None else np.full(len(PARAMS), np.nan)
    for j, p in enumerate(PARAMS):
        lo, hi = (to_phys(np.clip(u - sd_u * e, 0, 1))[j] for e in (1, -1)) if np.isfinite(u[j]) else (np.nan, np.nan)
        hi, lo = max(lo, hi), min(lo, hi)
        flag = "  <- determined by the prior" if data[j] < 0.5 else ""
        L.append(f"  {p.name:>11}: data share {data[j]:.2f}, ±1σ [{lo:.4g}, {hi:.4g}] (σ_u = {sd_u[j]:.2f}){flag}")
    return L


def save_jacobian(J):
    pd.DataFrame(J, index=[t.label for t in TARGETS], columns=[p.name for p in PARAMS]).to_csv(
        CAL_DIR / "jacobian.csv", float_format="%.4g")


# --------------------------------------------------------------------------
# Solving
# --------------------------------------------------------------------------
def calibrate(residual):
    u0 = np.clip(to_unit([spec_nominal(p) for p in PARAMS]), 0.02, 0.98)
    try:
        res = least_squares(residual, u0, jac=residual.jac, bounds=(0., 1.), method="trf",
                            x_scale=1.0, max_nfev=MAX_NFEV, verbose=1)
        F_end = res.fun[:len(TARGETS)]
        status = (f"finished: minimum of targets and priors found, {int(np.sum(np.abs(F_end) <= 1))} of {len(TARGETS)} "
                  f"targets within tolerance (least_squares: {res.message})")
        print(f"\nleast_squares: {res.message}")
        print("Not all targets within tolerance: noise too large, target at the bound of the search range or not "
              "reachable with these parameters -> see report.txt, jacobian.csv and evals.csv.")
    except Converged:
        status = "CONVERGED: all targets within tolerance"
        print("\nConverged: all targets within tolerance.")

    if residual.best is None:
        raise SystemExit("No point with all outputs defined, see CHECKS.")
    _, u, F = residual.best
    x = to_phys(u)
    write_params(x, CAL_DIR / "best_params.xml")
    print("\nBest point:", ", ".join(f"{p.name} = {xj:.4g}" for p, xj in zip(PARAMS, x)))
    print("  F =", ", ".join(f"{t.label} {Fi:+.2f}" for t, Fi in zip(TARGETS, F)))
    at_edge = [p.name for p, uj in zip(PARAMS, u) if uj <= 1e-6 or uj >= 1 - 1e-6]
    if at_edge:
        print("  At the bound of the search range:", ", ".join(at_edge), "-> widen the range?")
    print(f"  -> {CAL_DIR / 'best_params.xml'}")
    write_report(residual, status)
    return x


def check(procs):
    print(f"{len(PARAMS)} parameters, {len(TARGETS)} targets, {N_SEEDS} seeds, {HOURS} h after t*, template {config.TEMPLATE_XML.name}")
    print("Fixed at the nominal value (config.PARAMETERS, calibrate=False): "
          + (", ".join(f"{p.name} = {spec_nominal(p):.4g}" for p in FIXED_PARAMS) or "none"))
    for p in PARAMS:
        nom = spec_nominal(p)
        u0 = to_unit([spec_nominal(q) for q in PARAMS])[PARAMS.index(p)]
        sts = "+".join(config.SUBTYPE_NAMES[s] for s in p.subtypes)
        s = prior_spec(p)
        if s is None:
            pr = "no prior"
        else:
            c, w = (np.exp(s[0]), f"×/÷ {np.exp(s[1]):.3g}") if p.log else (s[0], f"± {s[1]:.3g}")
            pr = f"prior {s[2]}: {c:.4g} {w}" + (f" ({p.source})" if p.source else "")
        print(f"  {p.name:>11} ({sts:>10} {'='.join((p.param,) + p.also)}): [{p.low:.4g}, {p.high:.4g}]"
              f"{' log' if p.log else ''}, nominal {nom:.4g} (u = {u0:.2f}) | {pr}")
    for t in TARGETS:
        print(f"  target {t.label:>18} = {t.value:g} ± {t.tol_text}  ({t.source})")
    rounds = math.ceil(len(PARAMS) * N_SEEDS / procs)
    print(f"Per iteration: 1 point + {len(PARAMS)} difference points = {(len(PARAMS) + 1) * N_SEEDS} runs, "
          f"with {procs} cores 1 + {rounds} rounds of ~0.75 h")
    n_prior = sum(1 for p in PARAMS if prior_spec(p))
    print(f"Equations: {len(TARGETS)} targets + {n_prior} priors for {len(PARAMS)} parameters")
    if len(TARGETS) + n_prior < len(PARAMS):
        print(f"\nUNDERDETERMINED: the solution is not unique. Add targets or priors "
              f"(WEAK_PRIOR_SIGMA) or use --allow-underdetermined.")
        return False
    return True


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--allow-underdetermined", action="store_true")
    ap.add_argument("--seeds", type=int, default=N_SEEDS)
    ap.add_argument("--hours", type=float, default=HOURS)
    ap.add_argument("--max-nfev", type=int, default=MAX_NFEV)
    ap.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count())))
    ap.add_argument("--full", action="store_true", help="do not stop at the first point within tolerance, go on to the minimum")
    a = ap.parse_args()
    N_SEEDS, HOURS, MAX_NFEV = a.seeds, int(a.hours), a.max_nfev
    STOP_AT_TOLERANCE = not a.full
    determined = check(a.procs)
    if not a.check:
        if not determined and not a.allow_underdetermined:
            raise SystemExit(1)
        CAL_DIR.mkdir(exist_ok=True)
        with Pool(a.procs) as pool:
            r = Residual(pool, N_SEEDS, HOURS)
            def on_term(*_):   # SLURM sends SIGTERM at the time limit
                raise SystemExit("SIGTERM received, probably the job's time limit. Resubmit, finished runs are reused.")
            signal.signal(signal.SIGTERM, on_term)
            try:
                calibrate(r)
            except BaseException as e:
                msg = str(e) if isinstance(e, SystemExit) else traceback.format_exc()
                write_report(r, f"ABORTED: {msg}")
                raise
