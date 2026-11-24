"""
Comparison of two Morris analyses with the same design, e.g. old model (tree exclusion)
vs. new model (path length L = anaMinDist). Only reads finished results, computes nothing.

    sbatch sa_run.sbatch compare_morris.py                               # results_treeexcl vs. results
    sbatch sa_run.sbatch compare_morris.py results_treeexcl results      # old, new explicitly

Inputs per folder: morris_indices.csv (analyze_morris.py) and outputs_morris.csv (collect.py)
Outputs in results/ (or in the "new" folder):
  compare_morris.csv              mu*, mu*_rel, sigma, noise flag old/new per (time, output, factor)
  plots/morris_compare.png        heatmap mu*_rel old | new | new - old, dot = within noise
  plots/morris_compare_levels.png median of the outputs old vs. new per time
Console:
  1. per factor: in how many (output, time) above the noise, old vs. new, and largest mu*_rel
  2. per output and time: the three most important factors old and new
  3. output levels: median over all runs old, new and ratio (same factors and seeds)
"""
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import config

SA_OUTPUTS = config.SNAPSHOT_OUTPUTS + config.WAVE_OUTPUTS
INK, INK2, GRID = "#1f1f1e", "#5f5e5a", "#e4e3df"
SEQ = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
DIV = ["#b5541c", "#eb6834", "#f6b48f", "#fcfcfb", "#9ec5f4", "#3987e5", "#184f95"]  # new < old | new > old
OLD_C, NEW_C = "#8a8a85", "#2a78d6"


def load(folder):
    d = config.ROOT / folder
    ind = pd.read_csv(d / "morris_indices.csv")
    ind = ind[ind.region == "all"]
    out = pd.read_csv(d / "outputs_morris.csv")
    out = out[(out.region == "all") & out.output.isin(SA_OUTPUTS)]
    return ind, out


def columns(df):
    """(output, time) in the order of SA_OUTPUTS, existing ones only."""
    times = sorted(df.time.unique())
    have = set(zip(df.output, df.time))
    return [(o, t) for o in SA_OUTPUTS for t in times if (o, t) in have]


def plot_heatmaps(cmp, path):
    cols = columns(cmp)
    order = cmp.groupby("factor")[["rel_old", "rel_new"]].mean().max(axis=1).sort_values(ascending=False).index
    cmp = cmp.assign(col=list(zip(cmp.output, cmp.time)))
    piv = lambda v: cmp.pivot_table(index="factor", columns="col", values=v).reindex(index=order, columns=cols)
    panels = [("old (tree exclusion)", piv("rel_old"), piv("noise_old"), LinearSegmentedColormap.from_list("s", SEQ), 0, 1),
              ("new (path length L)", piv("rel_new"), piv("noise_new"), LinearSegmentedColormap.from_list("s", SEQ), 0, 1),
              ("new − old", piv("rel_new") - piv("rel_old"), None, LinearSegmentedColormap.from_list("d", DIV), -1, 1)]
    h = 0.3 * len(order) + 0.6
    fig, axes = plt.subplots(3, 1, figsize=(0.32 * len(cols) + 3, 3 * h + 2.2), sharex=True)
    for i, (ax, (title, M, N, cmap, lo, hi)) in enumerate(zip(axes, panels)):
        im = ax.imshow(M.values.astype(float), cmap=cmap, vmin=lo, vmax=hi, aspect="auto")
        if N is not None:   # dot: within noise
            yy, xx = np.where(N.values.astype(float) > 0.5)
            ax.scatter(xx, yy, s=6, color=INK2, lw=0)
        ax.set_xticks(range(len(cols)), [f"{o} @ {t:g} h" for o, t in cols], rotation=60, ha="right", fontsize=7, color=INK2)
        ax.tick_params(axis="x", labelbottom=(i == 2))
        ax.set_yticks(range(len(order)), order, fontsize=8, color=INK)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
        ax.set_title(title, loc="left", fontsize=10, color=INK)
        cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
        cb.outline.set_visible(False)
        cb.ax.tick_params(labelsize=7, colors=INK2)
    fig.suptitle("μ* / max μ* per output and time — dot: within noise (μ* ≤ 2 × noise level)",
                 x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_levels(lv, path):
    outs = [o for o in SA_OUTPUTS if o in set(lv.output)]
    n = len(outs); nc = min(n, 3); nr = int(np.ceil(n / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(4.2 * nc, 3.2 * nr + 0.5), squeeze=False)
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    for ax, o in zip(axes.flat, outs):
        d = lv[lv.output == o].sort_values("time")
        for lab, c, m, q1, q3 in (("old", OLD_C, "med_old", "q25_old", "q75_old"), ("new", NEW_C, "med_new", "q25_new", "q75_new")):
            ax.fill_between(d.time, d[q1], d[q3], color=c, alpha=0.15, lw=0)
            ax.plot(d.time, d[m], "o-", color=c, ms=4, lw=1.5, label=lab)
        ax.set_title(o, loc="left", fontsize=10, color=INK)
        ax.set_xlabel("h after t*", fontsize=8, color=INK2)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.grid(color=GRID, lw=0.6)
    axes.flat[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Output levels over all Morris runs: median and quartiles", x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    old_f, new_f = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ("results_treeexcl", "results")
    ind_old, out_old = load(old_f)
    ind_new, out_new = load(new_f)
    key = ["time", "output", "factor"]
    keep = ["mu_star", "mu_star_rel", "sigma", "in_noise", "n_valid"]
    cmp = ind_old[key + keep].merge(ind_new[key + keep], on=key, suffixes=("_old", "_new"))
    cmp = cmp.rename(columns={"mu_star_rel_old": "rel_old", "mu_star_rel_new": "rel_new",
                              "in_noise_old": "noise_old", "in_noise_new": "noise_new"})
    cmp["noise_old"] = cmp["noise_old"].astype(bool)
    cmp["noise_new"] = cmp["noise_new"].astype(bool)
    cmp["rel_diff"] = cmp["rel_new"] - cmp["rel_old"]
    res = config.ROOT / new_f
    (res / "plots").mkdir(exist_ok=True)
    cmp.to_csv(res / "compare_morris.csv", index=False)

    pd.set_option("display.width", 200)
    n_cells = cmp.groupby("factor").size().iloc[0]
    print(f"=== 1. Per factor: number of (output, time) above the noise (of {n_cells}) and largest μ*_rel ===")
    f = cmp.groupby("factor").agg(
        signif_old=("noise_old", lambda s: int((~s).sum())), signif_new=("noise_new", lambda s: int((~s).sum())),
        max_rel_old=("rel_old", "max"), max_rel_new=("rel_new", "max"),
        mean_rel_old=("rel_old", "mean"), mean_rel_new=("rel_new", "mean"))
    f["signif_diff"] = f.signif_new - f.signif_old
    print(f.sort_values("mean_rel_new", ascending=False).round(2).to_string())

    print("\n=== 2. Per output and time: most important factors (μ*_rel, * = within noise) ===")
    def top3(d, rel, noise):
        d = d.sort_values(rel, ascending=False).head(3)
        return ", ".join(f"{r.factor} {getattr(r, rel):.2f}{'*' if getattr(r, noise) else ''}" for r in d.itertuples())
    for o, t in columns(cmp):
        d = cmp[(cmp.output == o) & (cmp.time == t)]
        print(f"{o} @ {t:g} h")
        print(f"   old: {top3(d, 'rel_old', 'noise_old')}")
        print(f"   new: {top3(d, 'rel_new', 'noise_new')}")

    print("\n=== 3. Output levels: median over all runs (same design) ===")
    def stats(out, tag):
        g = out.groupby(["output", "time"])["value"]
        return pd.DataFrame({f"med_{tag}": g.median(), f"q25_{tag}": g.quantile(.25), f"q75_{tag}": g.quantile(.75)})
    lv = stats(out_old, "old").join(stats(out_new, "new"), how="inner").reset_index()
    lv["ratio"] = lv.med_new / lv.med_old.replace(0, np.nan)
    lv = lv.set_index(["output", "time"]).loc[[c for c in columns(lv)]].reset_index()
    print(lv[["output", "time", "med_old", "med_new", "ratio"]].round(4).to_string(index=False))
    # paired per run: ratio new/old with the same factors and seed
    p = out_old.merge(out_new, on=["run_id", "time", "output"], suffixes=("_old", "_new"))
    p = p[(p.value_old > 0) & (p.value_new > 0)]
    pr = p.assign(r=p.value_new / p.value_old).groupby(["output", "time"])["r"].median().rename("ratio_paired")
    print("\nMedian of the ratio new/old per run (paired):")
    print(pr.round(3).to_string())

    plot_heatmaps(cmp, res / "plots" / "morris_compare.png")
    plot_levels(lv, res / "plots" / "morris_compare_levels.png")
    print(f"\n-> {res / 'compare_morris.csv'}\n-> {res / 'plots' / 'morris_compare.png'}\n-> {res / 'plots' / 'morris_compare_levels.png'}")


if __name__ == "__main__":
    main()
