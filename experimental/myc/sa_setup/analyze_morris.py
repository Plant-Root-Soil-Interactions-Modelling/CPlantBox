"""
Phase 1: Morris analysis for all outputs and times (h after t*).

    python analyze_morris.py            # hyphal factors (design_morris.csv)
    python analyze_morris.py inoculum   # root/inoculum factors (design_inoculum.csv) -> results/inoculum_indices.csv,
                                        # plots/inoculum_*.png; additionally wait_days and the supply of hyphae (influx_*)
      -> results/morris_indices.csv     mu*, CI, mu, sigma per (time, region, output, factor)
      -> results/plots/morris_ranking_t<T>.png    mu* with 95% CI per output
      -> results/plots/morris_musigma_t<T>.png    mu* vs. sigma (nonlinearity/interaction)
      -> results/plots/morris_overview.png        mu*_rel for all outputs and times

mu* is given in output units per traversal of the whole factor range
(normalised factors) and is therefore comparable between factors, but not between outputs.
For comparison, mu*_rel = mu* / max(mu*) within one output is used.
mu*, mu, sigma and the bootstrap CI come from SALib.analyze.morris; the noise level
from the replicates is computed here.
"""
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from SALib.analyze import morris as morris_analyze

import config
from design import factors_of, morris_problem

KEY = ["time", "region", "output"]
# set by setup(tag)
TAG, FACTORS, NAMES, ROLE, GROUP, SA_OUTPUTS, PER_RUN, TRAJECTORIES, PREFIX = "morris", None, None, None, None, None, None, None, ""


def _group(f):
    if f.name == "hyphae_layer":
        return "geometry"
    if f.target != "hyphae":
        return "inoculum"
    return "control" if f.role == "control" else ("shared" if len(f.subtypes) > 1 else
                                                    config.SUBTYPE_NAMES[f.subtypes[0]].lower())


def setup(tag):
    global TAG, FACTORS, NAMES, ROLE, GROUP, SA_OUTPUTS, PER_RUN, TRAJECTORIES, PREFIX
    TAG, FACTORS = tag, factors_of(tag)
    NAMES = [f.name for f in FACTORS]
    ROLE = {f.name: f.role for f in FACTORS}
    GROUP = {f.name: _group(f) for f in FACTORS}
    # per run (replicate mean) and pooled per design point (calibration targets, evaluate_pooled.py)
    PER_RUN = list(dict.fromkeys(config.SNAPSHOT_OUTPUTS + config.IMAGE_OUTPUTS + config.WAVE_OUTPUTS + config.STATE_OUTPUTS
                                 + (config.INOCULUM_OUTPUTS if tag == "inoculum" else [])))
    SA_OUTPUTS = PER_RUN + config.POOLED_OUTPUTS
    TRAJECTORIES = config.MORRIS_TRAJECTORIES if tag == "morris" else config.INOCULUM_TRAJECTORIES
    PREFIX = "" if tag == "morris" else tag + "_"


setup("morris")

# colours (categories runner/BAS/shared/control, sequence for the heatmap); inoculum is never shown together with shared
COLORS = {"runner": "#2a78d6", "bas": "#eb6834", "shared": "#1b8a6b", "geometry": "#8a5fbf", "control": "#8a8a85",
          "inoculum": "#1b8a6b"}
LABELS = {"runner": "Runner", "bas": "BAS", "shared": "shared", "geometry": "geometry (layer)", "control": "control",
          "inoculum": "root / inoculum"}
MARKERS = {"runner": "o", "bas": "s", "shared": "D", "geometry": "P", "control": "X", "inoculum": "D"}
INK, INK2, GRID = "#1f1f1e", "#5f5e5a", "#e4e3df"
SEQ = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def _color(name):
    return COLORS[GROUP[name]]


def _grid(n):
    ncols = min(n, 3)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.4 * ncols, 4.2 * nrows + 0.6), squeeze=False)
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    return fig, axes.flat


def _elementary_effects(X, Y):
    """EEs per trajectory (only for the noise level, SALib does not return them).
    X: (r, k+1, k) normalised, Y: (..., r, k+1). Returns EE: (..., r, k), normalised as
    SALib: output change per step Δ in the normalised factor range."""
    dX = np.diff(X, axis=1)                              # (r, k, k): one factor changes per step
    i = np.abs(dX).argmax(-1)                            # (r, k): which factor in which step
    dx = np.take_along_axis(dX, i[..., None], -1)[..., 0]
    ee_step = np.diff(Y, axis=-1) / dx                   # (..., r, k) in step order
    EE = np.empty_like(ee_step)
    np.put_along_axis(EE, np.broadcast_to(i, ee_step.shape), ee_step, -1)
    return EE


NOISE_QUANTILE = 0.95
NOISE_DRAWS = 4000


def noise_threshold(se, seed=0):
    """Threshold for "above the noise" per factor (D21): μ* of a factor WITHOUT effect is the mean over the trajectories of
    |EE| with EE ~ N(0, se_t²) (se_t: standard error of the elementary effect of trajectory t, from the replicates or the
    bootstrap, so the uneven noise across the factor space is kept). Monte Carlo of this mean, NOISE_QUANTILE quantile.
    se: (trajectories, k) with nan for missing. Returns (mean level = E[μ*_noise], threshold)."""
    rng = np.random.default_rng(seed)
    k = se.shape[1]
    level, thr = np.full(k, np.nan), np.full(k, np.nan)
    for j in range(k):
        s = se[:, j][np.isfinite(se[:, j])]
        if len(s) == 0:
            continue
        draws = np.abs(rng.normal(0., 1., (NOISE_DRAWS, len(s))) * s).mean(1)
        level[j] = s.mean() * np.sqrt(2 / np.pi)
        thr[j] = np.quantile(draws, NOISE_QUANTILE)
    return level, thr


def prior_scale(factors):
    """Width of the prior (±2σ interval) relative to the Morris range, in the normalised (log) space; nan without a
    literature prior. μ* × this = effect of the actual uncertainty of the parameter instead of the whole range (D21)."""
    out = []
    for f in factors:
        p = getattr(f, "prior", None)
        if p is None:
            out.append(np.nan)
            continue
        g = np.log if f.log else (lambda v: v)
        out.append(abs(g(p[1]) - g(p[0])) / abs(g(f.high) - g(f.low)))
    return np.array(out)


def pooled_indices(problem, X):
    """Morris indices of the calibration targets pooled per design point (evaluate_pooled.py, pooled_<tag>.csv): one value
    per point, a trajectory counts if all k+1 points are defined. Noise level from the bootstrap SE of the points:
    SE of an EE = sqrt(se_s² + se_s+1²) / |Δ| (conservative, without the gain from shared seeds), × sqrt(2/π)."""
    path = config.RESULTS_DIR / f"pooled_{TAG}.csv"
    if not path.exists():
        print(f"No {path.name} (python evaluate_pooled.py results/design_{TAG}.csv): pooled calibration targets skipped.")
        return []
    k = problem["num_vars"]
    r, kp1, _ = X.shape
    dX = np.diff(X, axis=1)                                  # (r, k, k)
    fi = np.abs(dX).argmax(-1)                               # (r, k): factor changed in each step
    dx = np.abs(np.take_along_axis(dX, fi[..., None], -1)[..., 0])
    pl = pd.read_csv(path)
    rows = []
    for key, s in pl.groupby(KEY):
        Y, SE = np.full((r, kp1), np.nan), np.full((r, kp1), np.nan)
        Y[s.traj.values, s.step.values] = s.value.values
        SE[s.traj.values, s.step.values] = s.se.values
        valid = np.isfinite(Y).all(-1)
        n_valid = int(valid.sum())
        res = {m: np.full(k, np.nan) for m in ("mu_star", "mu", "sigma", "mu_star_lo", "mu_star_hi", "noise_level",
                                                "noise_threshold")}
        res["n_valid"] = np.full(k, n_valid)
        if n_valid:
            Si = morris_analyze.analyze(problem, X[valid].reshape(-1, k), Y[valid].reshape(-1), num_resamples=1000,
                                        conf_level=0.95, num_levels=config.MORRIS_LEVELS, seed=config.BASE_SEED)
            for m in ("mu_star", "mu", "sigma"):
                res[m] = np.asarray(Si[m], dtype=float)
            res["mu_star_lo"] = np.maximum(res["mu_star"] - Si["mu_star_conf"], 0)
            res["mu_star_hi"] = res["mu_star"] + Si["mu_star_conf"]
            se_ee = np.sqrt(SE[valid, 1:] ** 2 + SE[valid, :-1] ** 2) / dx[valid]     # (n_valid, k) in step order
            per_factor = np.full((n_valid, k), np.nan)
            np.put_along_axis(per_factor, fi[valid], se_ee, -1)
            res["noise_level"], res["noise_threshold"] = noise_threshold(per_factor)
        for j, n in enumerate(NAMES):
            rows.append(dict(zip(KEY, key), factor=n, role=ROLE[n], kind="pooled", **{m: float(res[m][j]) for m in res}))
    return rows


def compute():
    design = pd.read_csv(config.RESULTS_DIR / f"design_{TAG}.csv")
    problem = morris_problem(FACTORS)
    k = problem["num_vars"]
    X = np.load(config.RESULTS_DIR / f"{TAG}_X.npy").reshape(-1, k + 1, k)   # (r, k+1, k)
    out = pd.read_csv(config.RESULTS_DIR / f"outputs_{TAG}.csv").merge(
        design[["run_id", "traj", "step", "rep"]], on="run_id")
    out = out[out.output.isin(PER_RUN)]          # diagnostics ("n_*") and the per-run fits do not enter the SA
    r, kp1, _ = X.shape
    R = int(design["rep"].max()) + 1
    vals = out.set_index(KEY + ["rep", "traj", "step"])["value"]
    rows = []
    for key, s in vals.groupby(level=[0, 1, 2]):
        Y = np.full((R, r, kp1), np.nan)
        for (_, _, _, rep, t, st), v in s.items():
            Y[rep, t, st] = v
        # Every replicate is a self-consistent trajectory (shared seed) -> counts only when
        # complete. Valid are trajectories with at least one complete replicate.
        complete = np.isfinite(Y).all(-1)                  # (R, r)
        Y = np.where(complete[..., None], Y, np.nan)
        valid = complete.any(0)
        n_valid = int(valid.sum())
        res = {m: np.full(k, np.nan) for m in ("mu_star", "mu", "sigma")}
        res["n_valid"] = np.full(k, n_valid)
        res["mu_star_lo"], res["mu_star_hi"] = np.full(k, np.nan), np.full(k, np.nan)
        if n_valid:
            # EE is linear in Y: averaging replicates first == averaging EEs over replicates
            Y_mean = np.nanmean(Y[:, valid], 0)
            Si = morris_analyze.analyze(problem, X[valid].reshape(-1, k), Y_mean.reshape(-1),
                                        num_resamples=1000, conf_level=0.95,
                                        num_levels=config.MORRIS_LEVELS, seed=config.BASE_SEED)
            for m in ("mu_star", "mu", "sigma"):
                res[m] = np.asarray(Si[m], dtype=float)
            # SALib returns the CI as half-width (normal approximation of the bootstrap)
            res["mu_star_lo"] = np.maximum(res["mu_star"] - Si["mu_star_conf"], 0)
            res["mu_star_hi"] = res["mu_star"] + Si["mu_star_conf"]
        # Noise per factor: the mu* an ineffective factor would have from seed noise alone. SE of the replicate mean per
        # EE (valid trajectories); mean level E|N(0,se)| and the NOISE_QUANTILE threshold by Monte Carlo (noise_threshold)
        if R > 1:
            EE_rep = _elementary_effects(X, Y)             # (R, r, k)
            with np.errstate(invalid="ignore"):
                se = np.nanstd(EE_rep, 0, ddof=1) / np.sqrt(np.sum(np.isfinite(EE_rep), 0))
            res["noise_level"], res["noise_threshold"] = noise_threshold(se[valid])
        else:
            res["noise_level"], res["noise_threshold"] = np.full(k, np.nan), np.full(k, np.nan)
        for j, n in enumerate(NAMES):
            rows.append(dict(zip(KEY, key), factor=n, role=ROLE[n], kind="per run",
                             **{m: float(res[m][j]) for m in res}))
    rows += pooled_indices(problem, X)
    ind = pd.DataFrame(rows)
    mx = ind.groupby(KEY)["mu_star"].transform("max").replace(0, np.nan)
    ind["mu_star_rel"] = ind["mu_star"] / mx
    # "within noise": mu* not above the NOISE_QUANTILE quantile of what seed noise alone produces (Monte Carlo from the
    # per-trajectory standard errors, D21). Note: with shared seeds the control factor gives exactly 0 if it really has
    # no effect — it checks the b_prob switch, not the noise level.
    ind["in_noise"] = ~(ind["mu_star"] > ind["noise_threshold"])
    # effect of the actual uncertainty (prior width instead of the whole range) for parameters with a literature prior
    scale = dict(zip(NAMES, prior_scale(FACTORS)))
    ind["mu_star_prior"] = ind["mu_star"] * ind["factor"].map(scale)
    ind["mu_star_prior_rel"] = ind["mu_star_prior"] / mx      # relative to the largest μ* (whole range) of the output
    ind.to_csv(config.RESULTS_DIR / ("morris_indices.csv" if TAG == "morris" else f"{TAG}_indices.csv"), index=False)
    return ind


def _outputs_at(ind, t):
    present = set(ind[ind.time == t].output)
    return [o for o in SA_OUTPUTS if o in present]


def _legend_handles(kind):
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    groups = [g for g in LABELS if g in GROUP.values()]
    if kind == "bar":
        return [Patch(color=COLORS[g], label=LABELS[g]) for g in groups] + [
            Line2D([], [], marker="|", ls="", ms=10, mew=1.5, color=INK,
                   label=f"noise threshold ({NOISE_QUANTILE:.0%} quantile of μ* without effect)")]
    return [Line2D([], [], marker=MARKERS[g], ls="", color=COLORS[g], label=LABELS[g]) for g in groups] + [
        Line2D([], [], color=GRID, label="σ = μ*")]


def plot_ranking(ind, t):
    outs = _outputs_at(ind, t)
    sub = ind[(ind.time == t) & (ind.region == "all")]
    fig, axes = _grid(len(outs))
    for ax, o in zip(axes, outs):
        d = sub[sub.output == o].sort_values("mu_star")
        y = np.arange(len(d))
        ax.barh(y, d.mu_star, color=[_color(n) for n in d.factor], height=0.6, edgecolor="white", lw=1)
        ax.errorbar(d.mu_star, y, xerr=[d.mu_star - d.mu_star_lo, d.mu_star_hi - d.mu_star],
                    fmt="none", ecolor=INK2, elinewidth=0.8, capsize=2)
        if d.noise_threshold.notna().any():   # noise threshold per factor as a short tick
            ax.scatter(d.noise_threshold, y, marker="|", s=90, color=INK, lw=1.5, zorder=4)
        ax.set_yticks(y, d.factor, fontsize=8, color=INK)
        ax.set_title(o, fontsize=10, color=INK, loc="left")
        ax.set_xlabel("μ*", fontsize=8, color=INK2)
        _style(ax)
    fig.legend(handles=_legend_handles("bar"), loc="lower center", ncol=5, frameon=False, fontsize=9)
    fig.suptitle(f"Morris μ* with 95% bootstrap CI — t = {t:g} h after t*", x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    fig.savefig(config.RESULTS_DIR / "plots" / f"{PREFIX}morris_ranking_t{t:g}.png", dpi=150)
    plt.close(fig)


def plot_musigma(ind, t, n_label=6):
    outs = _outputs_at(ind, t)
    sub = ind[(ind.time == t) & (ind.region == "all")]
    fig, axes = _grid(len(outs))
    for ax, o in zip(axes, outs):
        d = sub[sub.output == o]
        for n, row in d.set_index("factor").iterrows():
            ax.scatter(row.mu_star, row.sigma, s=40, marker=MARKERS[GROUP[n]], color=_color(n),
                       edgecolor="white", lw=1, zorder=3)
        lim = np.nanmax([d.mu_star.max(), d.sigma.max(), 1e-12]) * 1.1
        ax.plot([0, lim], [0, lim], color=GRID, lw=1, zorder=1)        # σ = μ*
        for _, row in d.nlargest(n_label, "mu_star").iterrows():
            ax.annotate(row.factor, (row.mu_star, row.sigma), xytext=(4, 3),
                        textcoords="offset points", fontsize=7, color=INK)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_title(o, fontsize=10, color=INK, loc="left")
        ax.set_xlabel("μ*", fontsize=8, color=INK2); ax.set_ylabel("σ", fontsize=8, color=INK2)
        _style(ax); ax.grid(axis="y", color=GRID, lw=0.6)
    fig.legend(handles=_legend_handles("scatter"), loc="lower center", ncol=5, frameon=False, fontsize=9)
    fig.suptitle(f"Morris μ* vs. σ (large σ: nonlinear or interactions) — t = {t:g} h after t*",
                 x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    fig.savefig(config.RESULTS_DIR / "plots" / f"{PREFIX}morris_musigma_t{t:g}.png", dpi=150)
    plt.close(fig)


def plot_overview(ind):
    # all outputs and times in one heatmap: column = output @ time, row = factor
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("seq", SEQ)
    sub = ind[ind.region == "all"].copy()
    cols = [(o, t) for o in SA_OUTPUTS for t in sorted(sub.time.unique()) if ((sub.output == o) & (sub.time == t)).any()]
    sub["col"] = list(zip(sub.output, sub.time))
    order = sub.groupby("factor").mu_star_rel.mean().sort_values(ascending=False).index
    M = sub.pivot_table(index="factor", columns="col", values="mu_star_rel").reindex(index=order, columns=cols)
    fig, ax = plt.subplots(figsize=(0.55 * len(cols) + 2.5, 0.38 * len(order) + 2))
    im = ax.imshow(M.values, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(cols)), [f"{o} @ {t:g} h" for o, t in cols], rotation=60, ha="right", fontsize=8, color=INK2)
    ax.set_yticks(range(len(order)), order, fontsize=8, color=INK)
    for lbl in ax.get_yticklabels():
        lbl.set_color(_color(lbl.get_text()))
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("μ* / max μ* (per output and time)", fontsize=8, color=INK2)
    cb.outline.set_visible(False)
    ax.set_title("Relative sensitivity of all outputs (wave quantities over the whole window)", loc="left", fontsize=11, color=INK)
    fig.savefig(config.RESULTS_DIR / "plots" / f"{PREFIX}morris_overview.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


STATE_TARGETS = ["alpha", "beta", "K1", "K2", "front_velocity", "front_tip_speed", "bas_share_image", "L_RH_image",
                 "anastomoses", "root_efficiency_rel", "branch_angle_rh", "branch_angle_bas", "persistence_length"]


def state_correlation():
    """How strongly is each calibration target (per run, at max(SNAPSHOT_HOURS)) tied to the state of the colony at t*
    (D23)? Spearman rank correlation over all Morris runs with wait_days, n_crossings_tstar, influx_tstar and frozen_share.
    High |ρ|: the state at t* shapes the target; a parameter that acts strongly on that state (its μ* for the state
    outputs) then probably acts on the target indirectly. -> results/<tag>_state_correlation.csv"""
    path = config.RESULTS_DIR / f"outputs_{TAG}.csv"
    if not path.exists():
        return
    out = pd.read_csv(path)
    t_end = max(config.SNAPSHOT_HOURS)
    state = out[out.output.isin(config.STATE_OUTPUTS)].pivot_table(index="run_id", columns="output", values="value")
    targ = out[out.output.isin(STATE_TARGETS) & (out.time == t_end)].pivot_table(index="run_id", columns="output", values="value")
    if state.empty or targ.empty:
        return
    d = targ.join(state, how="inner")
    rows = []
    for t in [c for c in STATE_TARGETS if c in targ]:
        row = {"target": t}
        for s in [c for c in config.STATE_OUTPUTS if c in state]:
            x = d[[t, s]].dropna()
            row[s] = float(x[t].rank().corr(x[s].rank())) if len(x) > 10 else np.nan
        rows.append(row)
    tab = pd.DataFrame(rows).set_index("target")
    tab.to_csv(config.RESULTS_DIR / f"{TAG}_state_correlation.csv")
    print(f"\nCalibration targets (per run, {t_end:g} h) vs. the state of the colony at t*: Spearman ρ over all runs"
          "\n(|ρ| ≳ 0.5: the state at t* shapes the target; see the μ* of the parameters for wait_days, frozen_share, ...)")
    with pd.option_context("display.width", 160):
        print(tab.round(2).to_string())


CORR_FLAG = 0.9


def target_correlation():
    """E6: how much do the calibration targets say the same? Spearman rank correlation between the targets over the design
    points (mean of the replicates for the per-run targets, the pooled value for the pooled ones), each target at the
    times it is calibrated (calibrate.TARGETS). Strongly correlated targets count several times in the calibration:
    decide by hand whether one is enough or whether their tolerances are widened by √(number), as for the time courses.
    -> results/<tag>_target_correlation.csv, plots/<prefix>target_correlation.png"""
    try:
        import calibrate
    except Exception:
        return
    design = pd.read_csv(config.RESULTS_DIR / f"design_{TAG}.csv")[["run_id", "traj", "step"]]
    out_path, pooled_path = config.RESULTS_DIR / f"outputs_{TAG}.csv", config.RESULTS_DIR / f"pooled_{TAG}.csv"
    if not out_path.exists():
        return
    t_end = max(config.SNAPSHOT_HOURS)
    out = pd.read_csv(out_path).merge(design, on="run_id")
    pooled = pd.read_csv(pooled_path) if pooled_path.exists() else pd.DataFrame(columns=["traj", "step", "time", "output", "value"])
    cols = {}
    for t in calibrate.TARGETS:
        h = t.hour if t.hour is not None else t_end
        if t.output in config.POOLED_OUTPUTS and len(pooled):
            s = pooled[(pooled.output == t.output) & (pooled.time == t_end)].set_index(["traj", "step"]).value
        else:
            s = out[(out.output == t.output) & (out.time == h)].groupby(["traj", "step"]).value.mean()
        if len(s.dropna()) > 10:
            cols[t.label] = s
    if len(cols) < 2:
        return
    df = pd.DataFrame(cols)
    rho = df.rank().corr(min_periods=10)
    rho.to_csv(config.RESULTS_DIR / f"{TAG}_target_correlation.csv", float_format="%.3f")
    labels = list(rho.columns)
    pairs = [(a, b, rho.loc[a, b]) for i, a in enumerate(labels) for b in labels[i + 1:]
             if np.isfinite(rho.loc[a, b]) and abs(rho.loc[a, b]) >= CORR_FLAG]
    print(f"\nCalibration targets that say nearly the same (|Spearman ρ| ≥ {CORR_FLAG} over the design points):")
    if pairs:
        for a, b, r in sorted(pairs, key=lambda p: -abs(p[2])):
            print(f"  {a:26s} {b:26s} ρ = {r:+.2f}")
    else:
        print("  none")
    _plot_correlation(rho)


def _plot_correlation(rho):
    from matplotlib.colors import LinearSegmentedColormap
    # diverging: two hues, neutral grey in the middle (ρ = 0)
    cmap = LinearSegmentedColormap.from_list("div", [COLORS["runner"], "#f0efec", COLORS["bas"]])
    n = len(rho)
    fig, ax = plt.subplots(figsize=(0.42 * n + 2.5, 0.42 * n + 1.5))
    im = ax.imshow(rho.values, cmap=cmap, vmin=-1, vmax=1)
    ax.set_xticks(range(n), rho.columns, rotation=60, ha="right", fontsize=7, color=INK2)
    ax.set_yticks(range(n), rho.index, fontsize=7, color=INK)
    for i in range(n):
        for j in range(n):
            v = rho.values[i, j]
            if i != j and np.isfinite(v) and abs(v) >= CORR_FLAG:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6, color=INK)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("Spearman ρ over the design points", fontsize=8, color=INK2)
    cb.outline.set_visible(False)
    ax.set_title(f"Correlation of the calibration targets (values shown for |ρ| ≥ {CORR_FLAG})", loc="left", fontsize=11,
                 color=INK)
    fig.savefig(config.RESULTS_DIR / "plots" / f"{PREFIX}target_correlation.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def failure_map():
    """Where in the factor space do runs fail (D22)? Independent of SALib: per factor and Morris level (u = 0, 1/3, 2/3, 1)
    the share of the runs without result (no t*, aborted, missing) and the share of the design points whose pooled
    calibration targets are undefined. Failures that pile up at one end of a range bias μ* (those trajectories drop out).
    -> results/<tag>_failure_map.csv, plots/<prefix>morris_failure_map.png"""
    design = pd.read_csv(config.RESULTS_DIR / f"design_{TAG}.csv")
    status_path = config.RESULTS_DIR / f"status_{TAG}.csv"
    if not status_path.exists():
        return
    st = pd.read_csv(status_path)[["run_id", "status"]]
    d = design.merge(st, on="run_id", how="left")
    d["failed"] = d["status"].fillna("missing") != "ok"
    # excluded: no t* and ≤ EXCLUDE_CROSSED_CM hyphae on the fungal side, the exclusion criterion of the experiment
    d["excluded"] = d["status"] == "excluded"
    pooled_path = config.RESULTS_DIR / f"pooled_{TAG}.csv"
    undef = None
    if pooled_path.exists():   # a design point counts as undefined if any pooled target is nan there
        pl = pd.read_csv(pooled_path)
        undef = pl.groupby(["traj", "step"]).value.apply(lambda v: bool(v.isna().any())).rename("undefined").reset_index()
    rows = []
    levels = config.MORRIS_LEVELS - 1
    for n in NAMES:
        u = f"u_{n}"
        if u not in d:
            continue
        lv = np.round(d[u] * levels).astype(int)
        for L, g in d.groupby(lv):
            row = {"factor": n, "level": int(L), "u": L / levels, "n_runs": len(g), "failed_share": float(g.failed.mean()),
                   "excluded_share": float(g.excluded.mean())}
            if undef is not None:
                pts = g[["traj", "step"]].drop_duplicates().merge(undef, on=["traj", "step"], how="left")
                row["undefined_share"] = float(pts.undefined.fillna(True).mean())
            rows.append(row)
    tab = pd.DataFrame(rows)
    tab.to_csv(config.RESULTS_DIR / f"{TAG}_failure_map.csv", index=False)
    overall = float(d.failed.mean())
    print(f"\nFailures: {overall:.1%} of the runs without result, of them {d.excluded.mean():.1%} of all runs excluded as in "
          f"the experiment (no t*, ≤ {config.EXCLUDE_CROSSED_CM:g} cm hyphae on the fungal side)"
          + (f", design points with undefined pooled targets: {undef.undefined.mean():.1%}" if undef is not None else ""))
    spread = tab.groupby("factor").failed_share.agg(lambda s: s.max() - s.min()).sort_values(ascending=False)
    top = spread[spread > max(0.05, overall)]
    if len(top):
        print("Failures concentrated at one end of the range (max − min failed share over the levels):")
        for f, v in top.items():
            s = tab[tab.factor == f].set_index("level").failed_share
            print(f"  {f:12s} {v:.0%}  (levels 0..{levels}: " + ", ".join(f"{x:.0%}" for x in s) + ")")
    else:
        print("No factor with failures concentrated at one end of its range.")
    _plot_failure_map(tab, levels)


def _plot_failure_map(tab, levels):
    names = [n for n in NAMES if n in set(tab.factor)]
    if not names:
        return
    ncols = min(len(names), 6)
    nrows = int(np.ceil(len(names) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.2 * ncols, 2.0 * nrows + 0.8), squeeze=False, sharey=True)
    ymax = max(0.05, float(np.nanmax(tab[[c for c in ("failed_share", "undefined_share") if c in tab]].values)) * 1.1)
    for ax, n in zip(axes.flat, names):
        s = tab[tab.factor == n].sort_values("level")
        x = s.level.values
        ax.bar(x - 0.18, s.failed_share, width=0.34, color=COLORS[GROUP[n]], label="runs without result")
        ax.bar(x - 0.18, s.excluded_share, width=0.34, color=INK, alpha=0.55,
               label=f"of them excluded (≤ {config.EXCLUDE_CROSSED_CM:g} cm, as experiment)")
        if "undefined_share" in s:
            ax.bar(x + 0.18, s.undefined_share, width=0.34, color=INK2, alpha=0.45, label="points with undefined targets")
        ax.set_xticks(range(levels + 1), [f"{i / levels:.2g}" for i in range(levels + 1)], fontsize=7)
        ax.set_ylim(0, ymax)
        ax.set_title(n, fontsize=9, color=INK, loc="left")
        _style(ax)
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
    for ax in list(axes.flat)[len(names):]:
        ax.set_visible(False)
    h, l = axes.flat[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, frameon=False, fontsize=9)
    fig.suptitle("Failures per factor and Morris level (normalised position in the range)", x=0.01, ha="left",
                 fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    fig.savefig(config.RESULTS_DIR / "plots" / f"{PREFIX}morris_failure_map.png", dpi=150)
    plt.close(fig)


def main():
    (config.RESULTS_DIR / "plots").mkdir(parents=True, exist_ok=True)
    ind = compute()
    for t in sorted(ind.time.unique()):
        plot_ranking(ind, t); plot_musigma(ind, t)
    plot_overview(ind)
    failure_map()
    state_correlation()
    target_correlation()
    sub = ind[ind.region == "all"]
    top = (sub[~sub.in_noise].groupby("factor").mu_star_rel.max().sort_values(ascending=False))
    print("\nFactors above the noise level in at least one output and time,"
          "\nvalue = largest relative μ* over outputs and times:")
    print(top.round(2).to_string())
    rest = sorted(set(NAMES) - set(top.index))
    print(f"\nAbove the noise level in no output: {', '.join(rest) or '-'}")
    ctrl = sub[sub.role == "control"]
    if len(ctrl) and ctrl.mu_star.max() > 0:
        print("\nWarning: the control factor has μ* > 0 -> does the b_prob switch work as expected?")
    nv = ind.n_valid.min()
    if nv < TRAJECTORIES:
        print(f"\nWarning: because of failed runs, some outputs have only {nv} valid trajectories.")
    if TAG == "inoculum":
        compare_to_hyphae(ind)
    else:
        calibration_suggestion(ind)


SUGGEST_REL = 0.1   # suggestion "free": μ* above the noise level and at least this share of the largest μ* of a target


def calibration_suggestion(ind):
    """Suggestion (ONLY a suggestion, nothing is changed) which parameters to calibrate (config.PARAMETERS, calibrate=True):
    those whose μ* lies above the noise level and is at least SUGGEST_REL of the largest μ* for at least one calibration
    target (the pooled outputs and the image graph quantities that are targets, at max(SNAPSHOT_HOURS)). Set calibrate by
    hand in config.py."""
    t = max(config.SNAPSHOT_HOURS)
    targets = config.POOLED_OUTPUTS + ["root_efficiency_rel", "bas_share_image", "L_RH_image", "anastomoses"]
    sub = ind[(ind.region == "all") & (ind.time == t) & ind.output.isin(targets)]
    if sub.empty:
        return
    rows = []
    for p in config.PARAMETERS:
        s = sub[sub.factor == p.name]
        if s.empty:
            continue
        clear = s[(~s.in_noise) & (s.mu_star_rel >= SUGGEST_REL)]
        # with a literature prior: does its ACTUAL uncertainty (prior width) still matter? (D21)
        prior_matters = p.prior is None or bool((clear.mu_star_prior_rel >= SUGGEST_REL).any())
        suggestion = "free" if len(clear) and prior_matters else \
            ("fix (prior narrow: its uncertainty hardly matters)" if len(clear) else "fix at nominal")
        rows.append({"parameter": p.name, "calibrate_in_config": p.calibrate, "suggestion": suggestion,
                     "max_mu_star_rel": round(float(s.mu_star_rel.max()), 3),
                     "max_mu_star_prior_rel": round(float(s.mu_star_prior_rel.max()), 3) if p.prior else np.nan,
                     "targets": ", ".join(f"{o} ({r:.2f})" for o, r in clear.sort_values("mu_star_rel", ascending=False)
                                          [["output", "mu_star_rel"]].itertuples(index=False))})
    tab = pd.DataFrame(rows)
    tab.to_csv(config.RESULTS_DIR / "morris_calibration_suggestion.csv", index=False)
    print(f"\nSuggestion for the calibration (μ* above the noise threshold and ≥ {SUGGEST_REL:.0%} of the largest μ* for a "
          f"target at {t:g} h;\nwith a literature prior also the μ* scaled to the prior width ≥ {SUGGEST_REL:.0%});"
          "\nONLY a suggestion: set calibrate=True/False by hand in config.PARAMETERS. The parameters marked free are also"
          "\nthe candidates for the variance-based (Sobol) step after the calibration")
    with pd.option_context("display.width", 220, "display.max_colwidth", 90):
        print(tab.to_string(index=False))
    anastomosis_group(sub)


ANASTOMOSIS_GROUP = ["distTH", "ana", "ana_L_scale", "R_a", "B_a"]


def anastomosis_group(sub):
    """These parameters act on the anastomosis rate nearly interchangeably (the targets see mostly their product).
    μ* relative to distTH for the anastomosis targets: if distTH clearly dominates, fix ana_L_scale and the radii a (and
    possibly ana) at their nominal values in config.PARAMETERS (calibrate=False), so that one parameter carries the role."""
    s = sub[sub.output.isin(["alpha", "anastomoses"]) & sub.factor.isin(ANASTOMOSIS_GROUP)]
    if s.empty:
        return
    w = s.pivot_table(index="factor", columns="output", values="mu_star").reindex(ANASTOMOSIS_GROUP)
    ref = w.loc["distTH"] if "distTH" in w.index else None
    print("\nAnastomosis parameters (nearly interchangeable): μ* for the anastomosis targets"
          + (", relative to distTH" if ref is not None else ""))
    with pd.option_context("display.width", 160):
        print((w / ref if ref is not None else w).round(2).to_string())


def compare_to_hyphae(ind):
    """Decision for the root/inoculum factors: largest μ* of these factors relative to the largest μ* of the hyphal factors
    (morris_indices.csv), per output and time. μ* is in output units per traversal of the factor range in both, so the
    ratio says how much the outputs depend on the root phase compared with the hyphal parameters (within their ranges)."""
    path = config.RESULTS_DIR / "morris_indices.csv"
    if not path.exists():
        print("\nNo morris_indices.csv of the hyphal factors yet: comparison skipped.")
        return
    hy = pd.read_csv(path)
    hy = hy[(hy.region == "all") & (hy.role != "control")].groupby(["time", "output"]).mu_star.max().rename("mu_hyphae")
    ino = ind[ind.region == "all"].groupby(["time", "output"]).agg(mu_inoculum=("mu_star", "max"),
                                                                     top=("mu_star", lambda s: ind.loc[s.idxmax(), "factor"] if s.max() > 0 else ""))
    cmp = ino.join(hy, how="inner")
    cmp["ratio"] = cmp.mu_inoculum / cmp.mu_hyphae.replace(0, np.nan)
    cmp.reset_index().to_csv(config.RESULTS_DIR / "inoculum_vs_hyphae.csv", index=False)
    print("\nRoot/inoculum vs. hyphal factors: max μ* (inoculum) / max μ* (hyphae) per output "
          "(<< 1: the normalisation to t* removes the root phase; ~1 or more: root phase matters):")
    with pd.option_context("display.width", 200):
        print(cmp.round(3).to_string())


if __name__ == "__main__":
    setup(sys.argv[1] if len(sys.argv) > 1 else "morris")
    main()
