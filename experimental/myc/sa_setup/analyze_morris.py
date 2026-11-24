"""
Phase 1: Morris analysis for all outputs and times (h after t*).

    python analyze_morris.py
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
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from SALib.analyze import morris as morris_analyze

import config
from design import morris_problem

KEY = ["time", "region", "output"]
NAMES = [f.name for f in config.FACTORS]
ROLE = {f.name: f.role for f in config.FACTORS}
GROUP = {f.name: "control" if f.role == "control" else ("shared" if len(f.subtypes) > 1 else
         config.SUBTYPE_NAMES[f.subtypes[0]].lower()) for f in config.FACTORS}
SA_OUTPUTS = config.SNAPSHOT_OUTPUTS + config.WAVE_OUTPUTS

# colours (categories runner/BAS/shared/control, sequence for the heatmap)
COLORS = {"runner": "#2a78d6", "bas": "#eb6834", "shared": "#1b8a6b", "control": "#8a8a85"}
LABELS = {"runner": "Runner", "bas": "BAS", "shared": "shared", "control": "control"}
MARKERS = {"runner": "o", "bas": "s", "shared": "D", "control": "X"}
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


def compute():
    design = pd.read_csv(config.RESULTS_DIR / "design_morris.csv")
    problem = morris_problem()
    k = problem["num_vars"]
    X = np.load(config.RESULTS_DIR / "morris_X.npy").reshape(-1, k + 1, k)   # (r, k+1, k)
    out = pd.read_csv(config.RESULTS_DIR / "outputs_morris.csv").merge(
        design[["run_id", "traj", "step", "rep"]], on="run_id")
    out = out[out.output.isin(SA_OUTPUTS)]       # diagnostics ("n_*", wait_days) do not enter the SA
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
        # Noise level per factor: the mu* an ineffective factor would have from seed noise alone.
        # SE of the replicate mean per EE, E|N(0,se)| = se*sqrt(2/pi), averaged over trajectories.
        if R > 1:
            EE_rep = _elementary_effects(X, Y)             # (R, r, k)
            se = np.nanstd(EE_rep, 0, ddof=1) / np.sqrt(R)
            res["noise_level"] = np.nanmean(se, 0) * np.sqrt(2 / np.pi)
        else:
            res["noise_level"] = np.full(k, np.nan)
        for j, n in enumerate(NAMES):
            rows.append(dict(zip(KEY, key), factor=n, role=ROLE[n],
                             **{m: float(res[m][j]) for m in res}))
    ind = pd.DataFrame(rows)
    mx = ind.groupby(KEY)["mu_star"].transform("max").replace(0, np.nan)
    ind["mu_star_rel"] = ind["mu_star"] / mx
    # "within noise": mu* not clearly (factor 2) above what seed noise alone produces.
    # Note: with shared seeds the control factor gives exactly 0 if it really has no
    # effect — it checks the b_prob switch, not the noise level.
    ind["in_noise"] = ind["mu_star"] <= 2 * ind["noise_level"]
    ind.to_csv(config.RESULTS_DIR / "morris_indices.csv", index=False)
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
            Line2D([], [], marker="|", ls="", ms=10, mew=1.5, color=INK, label="2 × noise level (from replicates)")]
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
        if d.noise_level.notna().any():   # noise level per factor as a short tick
            ax.scatter(2 * d.noise_level, y, marker="|", s=90, color=INK, lw=1.5, zorder=4)
        ax.set_yticks(y, d.factor, fontsize=8, color=INK)
        ax.set_title(o, fontsize=10, color=INK, loc="left")
        ax.set_xlabel("μ*", fontsize=8, color=INK2)
        _style(ax)
    fig.legend(handles=_legend_handles("bar"), loc="lower center", ncol=5, frameon=False, fontsize=9)
    fig.suptitle(f"Morris μ* with 95% bootstrap CI — t = {t:g} h after t*", x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    fig.savefig(config.RESULTS_DIR / "plots" / f"morris_ranking_t{t:g}.png", dpi=150)
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
    fig.savefig(config.RESULTS_DIR / "plots" / f"morris_musigma_t{t:g}.png", dpi=150)
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
    fig.savefig(config.RESULTS_DIR / "plots" / "morris_overview.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    (config.RESULTS_DIR / "plots").mkdir(parents=True, exist_ok=True)
    ind = compute()
    for t in sorted(ind.time.unique()):
        plot_ranking(ind, t); plot_musigma(ind, t)
    plot_overview(ind)
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
    if nv < config.MORRIS_TRAJECTORIES:
        print(f"\nWarning: because of failed runs, some outputs have only {nv} valid trajectories.")


if __name__ == "__main__":
    main()
