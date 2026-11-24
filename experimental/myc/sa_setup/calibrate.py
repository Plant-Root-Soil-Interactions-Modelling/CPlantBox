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
from xml_tools import hyphae_nominal, write_hyphae_xml

amfd = ev.amfd
RUNNER, BAS = config.RUNNER, config.BAS
deg = math.radians


# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------
@dataclass
class CalParam:
    name: str
    subtypes: tuple             # subTypes of the <hyphae> the value is written to
    param: str                  # name of the <parameter>; "successor" = share of the first successor (runner: BAS)
    low: float                  # search range
    high: float
    log: bool = False
    also: tuple = ()            # further XML parameters of the same subTypes that get the same value (lb = ln)
    prior: tuple | None = None  # (lo, hi): literature range as ±2σ interval (log parameters: in log); None = weak prior
    source: str = ""            # origin of the prior


# Separate for runner/BAS as in the SA (config.FACTORS), additionally separate tropismN and a (radius: they look
# very different in the experiment, acts on anastomosis via the tip–hypha distance); otherwise shared.
# Not included: b of the runners (no effect with b_prob = 0), ln and lmax of the BAS (no laterals, lmax ≈ 3 dx),
# la (fixed, very small, template XML). lb of the runners is coupled to ln; lb of the BAS stays fixed.
# Ranges as in the SA where available; R_ln, B_b, B_hlt, distTH wider, tropismN, a, successor estimated.
# Priors from the sources in config.FACTORS. B_b without literature prior: the measurement (Oyarte Galvez 2025) is
# already the target branching_rate, it would count twice otherwise.
V_SRC = "180–300 µm/h around 240 µm/h (Oyarte Galvez et al. 2025)"
THETA_SRC = "60° (Schnepf 2016) to ≈ 70° (Oyarte Galvez 2025), widened to 55–75° (different angle definitions)"
PARAMS = [
    # ---- RunnerHyphae (subType 1) ----
    CalParam("R_v", (RUNNER,), "v", 0.432, 0.72, prior=(0.432, 0.72), source=V_SRC),
    CalParam("R_hlt", (RUNNER,), "hlt", 4.0, 10.0, prior=(4.0, 12.5),
             source="SA range Schnepf 2016 4–10 d, death rates Schnepf 2008 4–12.5 d"),
    CalParam("R_lmax", (RUNNER,), "lmax", 1.5, 5.0),
    CalParam("R_ln", (RUNNER,), "ln", 0.05, 0.6, log=True, also=("lb",), prior=(0.1, 0.6),
             source="~1 mm (Giovannetti, in Schnepf 2008) to ~6 mm (Oyarte Galvez 2025)"),
    CalParam("R_theta", (RUNNER,), "theta", deg(50), deg(80), prior=(deg(55), deg(75)), source=THETA_SRC),
    CalParam("R_tropismN", (RUNNER,), "tropismN", 1.0, 8.0),
    CalParam("R_succ_BAS", (RUNNER,), "successor", 0.5, 0.95),   # share of BAS among the runner laterals
    # radius [cm], nominal 0.005 = 50 µm (Schnepf et al. 2017, radius); lower bound of the search range 10 µm (Bago et al. 1998:
    # runner diameter > 20 µm). Prior 25–100 µm (factor 2 around 50 µm, log) as ±2σ: Schnepf gives only a point value
    CalParam("R_a", (RUNNER,), "a", 0.001, 0.01, log=True, prior=(0.0025, 0.01),
             source="radius 50 µm (Schnepf et al. 2017), ±2σ = factor 2"),
    # ---- BAS (subType 2) ----
    CalParam("B_v", (BAS,), "v", 0.432, 0.72, prior=(0.432, 0.72), source=V_SRC),
    CalParam("B_b", (BAS,), "b", 0.3, 2.0, log=True),
    CalParam("B_hlt", (BAS,), "hlt", 3.0, 12.0, prior=(5.0, 9.0), source="BAS degenerate after ~7 d (Bago 1998; Kokkoris 2026)"),
    CalParam("B_theta", (BAS,), "theta", deg(50), deg(80), prior=(deg(55), deg(75)), source=THETA_SRC),
    CalParam("B_tropismN", (BAS,), "tropismN", 0.5, 4.0),
    CalParam("B_a", (BAS,), "a", 0.00005, 0.005, log=True, prior=(0.0001, 0.0005),   # radius [cm], nominal 0.0025
             source="diameter 2–10 µm = radius 1–5 µm (Bago et al. 1998, Mycologia, doi 10.1080/00275514.1998.12026878)"),
    # ---- shared ----
    CalParam("ana", (RUNNER, BAS), "ana", 0.5, 1.0),
    CalParam("distTH", (RUNNER, BAS), "distTH", 0.005, 0.03, log=True),
]
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

    def residual(self, y):
        if not np.isfinite(y) or (self.log and y <= 0):
            return np.nan
        return np.log(y / self.value) / np.log(1 + self.tol_rel) if self.log else (y - self.value) / (self.tol_rel * self.value)

    @property
    def label(self):
        return self.output if self.hour is None else f"{self.output}@{self.hour:g}h"


# Available outputs. Wave quantities over the window 0..HOURS, pooled over all seeds:
WAVE_OUTPUTS = {
    "beta": "[µm/h] a = beta rho n with a, n of the runner tips (ratio of sums)",
    "branching_rate": "[h^-1] new runner tips per growing runner tip and hour",
    "encounter_C": "[mm^2 µm^-1 h^-1] f = C rho",
    "K1": "[µm mm^-2] saturation of the filament density (median of saturated rings)",
    "K2": "[mm^-2] peak of the runner tip density (median of rings with the peak in the window)",
    "wave_width": "[h] duration of the runner tip wave in a ring (median)",
    "front_velocity": "[mm/d] ring radius vs. arrival t_n of the filament density (sigmoid fit, Theil-Sen, all seeds)",
    "front_speed": "[mm/d] slope of the tip-weighted radius from 10 h (mean over seeds)",
    "root_efficiency_rel": "(E_r - E_r,MST) / (E_r,DT - E_r,MST), 2D, at the end of the simulation (mean over seeds)",
}
# Snapshot outputs (mean over seeds) at time Target.hour: config.SNAPSHOT_OUTPUTS
# total_length [cm], active_tips, bas_share, anastomoses (since t*), front_position [mm], spread_r90 [mm],
# max_tip_density [mm^-2]

# A unique solution needs at least as many targets as parameters (here 16), and every
# parameter must clearly affect at least one target. Add further measurements as Target(...).
TARGETS = [
    Target("beta", 23.0, 0.10),                    # Oyarte Galvez et al. 2025
    Target("K1", 1400.0, 0.10),                    # saturation density, Fig. 2b, dashed line
    Target("bas_share", 0.30, 0.10, hour=60),      # BAS share of the hyphal length
    Target("K2", 0.15, 0.33),                      # Fig. 2b, tip peak 0.1–0.2
    Target("branching_rate", 0.04, 0.20),          # 4 % per hour
    Target("front_velocity", 280e-3 * 24, 0.10),   # wave speed 280 µm/h = 6.72 mm/d
    Target("root_efficiency_rel", 0.5, 0.20, log=False),   # relative root efficiency, end (HOURS), can be < 0
    # Target("wave_width", ..., 0.20),
    # Target("total_length", ..., 0.15, hour=30),
    # Target("active_tips", ..., 0.15, hour=30),
    # Target("anastomoses", ..., 0.20, hour=60),
]

HOURS = config.HOURS_HYPHAE   # simulation time after t* [h], window of the wave quantities
N_SEEDS = 16                  # seeds per parameter point
SEED0 = config.BASE_SEED + 300_000
DIFF_STEP = 0.15              # difference step in the normalised space
MAX_NFEV = 30                 # max. solver iterations (least_squares: max_nfev, without differences)
MIN_OK = 0.5                  # share of seeds that must reach t*, otherwise the point counts as invalid
NAN_PENALTY = 10.             # F_i for undefined outputs (fit failed), so that the solver keeps going
N_BOOT = 200
CAL_DIR = config.ROOT / "calibration"
# Version of the raw data (raw.npz) in the key of the runs: increase when model_hook stores new quantities that the
# targets need; old runs are then not reused. 2: branch_new_runner
RAW_VERSION = 2


# --------------------------------------------------------------------------
# Parameter space
# --------------------------------------------------------------------------
def spec_nominal(p):
    return hyphae_nominal(p.subtypes[0], p.param)


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
    """Physical parameter values -> {(subType, name): value} for write_hyphae_xml."""
    upd = {}
    for p, xj in zip(PARAMS, x):
        for st in p.subtypes:
            for name in (p.param,) + p.also:
                upd[(st, name)] = xj
    return upd


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
    xml_path = write_hyphae_xml(updates, rdir / "params.xml")
    meta = {"seed": seed, "hours": hours, "values": {f"{st}:{n}": v for (st, n), v in updates.items()}}
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
    sim_rh = dict(sim, tip_dens=sim["active_tips_runner"] / np.asarray(sim["area"]).reshape(-1, 1))
    return {"sim": sim, "table": ev.fit_table(sim_rh), "front_speed": ev.front_speed(sim, HOURS),
            "snap": {h: ev.snapshot_outputs(sim, h) for h in hours}}


def pooled_outputs(preps):
    sims = [p["sim"] for p in preps]
    good, good_K1, good_tips = ev.good_rings(pd.concat([p["table"] for p in preps], ignore_index=True))
    n_min = config.MIN_RINGS
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wave = {
            # a = beta rho n with a and n of the runner tips (growing tips as in the experiment), rho of all hyphae
            "beta": amfd.anastomosisBeta(sims, min_tips=config.ENCOUNTER_MIN_TIPS)["runner"],
            # new runner tips per growing runner tip and hour (alpha in b(n) = alpha n)
            "branching_rate": amfd.branchingRate(sims, min_tips=config.ENCOUNTER_MIN_TIPS, suffix="_runner"),
            "encounter_C": amfd.encounterC(sims, min_tips=config.ENCOUNTER_MIN_TIPS)["total"],
            "K1": np.median(good_K1.K1) if len(good_K1) >= n_min else np.nan,
            "K2": np.median(good_tips.K2) if len(good_tips) >= n_min else np.nan,   # peak of the runner tip density
            "wave_width": np.median(ev.FWHM_BELL / good_tips.lam_tip) if len(good_tips) >= n_min else np.nan,
            # wave position from the sigmoid fit of the filament density (Fig. 2d inset): ring radius vs. arrival time t_n
            "front_velocity": 24 * amfd.theilSlope(good.t_n, good.r)[0] if len(good) >= n_min else np.nan,
            "front_speed": np.nanmean([p["front_speed"] for p in preps]),
            "root_efficiency_rel": np.nanmean([float(p["sim"].get("root_eff_E_rel", np.nan)) for p in preps]),
        }
    out = {}
    for t in TARGETS:
        if t.hour is None:
            out[t.label] = float(wave[t.output])
        else:
            out[t.label] = float(np.nanmean([p["snap"][t.hour][t.output] for p in preps]))
    return out


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
               "norm_F": float(np.linalg.norm(F)), "norm_prior": float(np.linalg.norm(z)), "cost": cost}
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
        if np.all(np.abs(F) <= 1):
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
        print(f"   {n:>18}: {rec[n]:10.4g} ± {rec[f'se_{n}']:.2g}  (target {t.value:g} ± {t.tol_rel:.0%})  "
              f"F = {rec[f'F_{n}']:+.2f}  {flag}")
    print(f"   |F| = {rec['norm_F']:.3f}, |Prior| = {rec['norm_prior']:.3f}", flush=True)


HINTS = {
    "K1": "too few saturated rings with good fits -> longer simulation (--hours 120) or more seeds",
    "K2": "too few rings with a complete tip peak -> longer simulation (--hours 120) or more seeds",
    "wave_width": "as K2", "front_velocity": "as K2",
    "root_efficiency_rel": "root efficiency not computed -> see root efficiency errors",
    "beta": "too few rings with >= ENCOUNTER_MIN_TIPS tips", "branching_rate": "as beta",
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
          f"Evaluations {len(recs)} ({n_pts} iteration points, {len(recs) - n_pts} difference points), runs {r.n_runs}", ""]

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
            L.append(f"  {t.label:>20}: {y:10.4g} ± {se:<8.2g} target {t.value:g} ± {t.tol_rel:.0%}   F = {Fi:+6.2f}  {flag}")
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
    if not any(l.startswith(("Failed", "Undef", "Noise", "Root", "No clear", "Conflict")) and not l.startswith("Failed runs: 0")
               for l in L[L.index("CHECKS"):]):
        L.append("No issues found.")
    (CAL_DIR / "report.txt").write_text("\n".join(L) + "\n", encoding="utf-8")


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
    write_hyphae_xml(xml_updates(x), CAL_DIR / "best_params.xml")
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
        print(f"  target {t.label:>18} = {t.value:g} ± {t.tol_rel:.0%}")
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
    a = ap.parse_args()
    N_SEEDS, HOURS, MAX_NFEV = a.seeds, int(a.hours), a.max_nfev
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
