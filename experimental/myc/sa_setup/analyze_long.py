"""
Test with a longer simulation (design.py long), at the nominal point over LONG_HOURS h:

1. Are the wave quantities converged at HOURS_HYPHAE = 120 h (compare 150, 180 h), are the wave fits defined in almost every run?
   For Morris an output must exist in all k+1 points of a trajectory, i.e. in almost all runs.
2. Roots turned off at t* (group "long", as in SA and calibration) vs. roots growing on (group "long_roots_on"),
   same seeds: when do new hyphae stop coming with the roots turned off, and from when on do the outputs differ?

    python analyze_long.py
      -> results/long_summary.csv          valid runs, median, CV per output, time window and group
      -> results/long_roots_compare.csv    roots on vs. off: relative difference per output and time, paired over seeds
      -> results/long_influx.csv           supply of hyphae per hour (mean over seeds per group)
      -> results/plots/long_influx.png
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

import config
import model_hook

SHOW = ["front_speed", "spread_r90", "max_tip_density", "front_velocity", "v_tips", "v_length_sat", "K1", "K2", "wave_width", "encounter_C",
        "alpha", "beta", "root_efficiency_rel",
        "n_wave_rings", "n_K1_rings", "n_tip_rings", "front_position", "total_length", "active_tips", "bas_share", "anastomoses",
        "wait_days"]   # wait_days: up to t* both groups are identical (same seeds), the difference must be 0
OFF, ON = "long", "long_roots_on"
LABELS = {OFF: "roots off at t*", ON: "roots growing on"}
COLORS = {OFF: "#2a78d6", ON: "#eb6834"}
INK, INK2, GRID = "#1f1f1e", "#5f5e5a", "#e4e3df"
# quantities of model_hook.influx_diagnostics: key in raw.npz -> label
INFLUX = {"diag_cross_in": "hyphae crossing into the fungal side [h⁻¹]",
          "diag_cross_in_runner": "of them runner [h⁻¹]",
          "diag_root_hyphae_new": "new hyphae from roots [h⁻¹]",
          "diag_root_length_new": "new root length [cm h⁻¹]"}
DROP = 0.1   # "supply stopped": 6 h mean of the crossings below this share of their maximum


def summary(out, status, design):
    k1 = len(config.FACTORS) + 1
    rows = []
    for g, d in design.groupby("group"):
        st = status[status.run_id.isin(d.run_id)]
        n = len(st)
        ok = int((st.status == "ok").sum())
        t = st.sim_seconds.dropna() / 3600
        print(f"{LABELS.get(g, g)}: {ok} of {n} runs ok, simulation median {t.median():.2f} h, max {t.max():.2f} h")
        sub = out[out.run_id.isin(d.run_id) & out.output.isin(SHOW)]
        for (o, h), s in sub.groupby(["output", "time"]):
            v = s.value
            p = v.notna().sum() / n if n else np.nan
            rows.append({"group": g, "output": o, "time": h, "valid": f"{v.notna().sum()}/{n}", "median": v.median(),
                         "cv": v.std() / v.mean() if v.mean() else np.nan,
                         # share of complete Morris replicates if every point has a value with probability p
                         "p_traj": p ** k1})
    tab = pd.DataFrame(rows)
    tab.to_csv(config.RESULTS_DIR / "long_summary.csv", index=False)
    off = tab[tab.group == OFF]
    if len(off):
        wide = off.pivot(index="output", columns="time", values=["valid", "median", "cv", "p_traj"])
        wide = wide.reindex([o for o in SHOW if o in wide.index])
        with pd.option_context("display.width", 220, "display.max_columns", 50):
            print(f"\n{LABELS[OFF]}, per output and time window (h after t*): valid runs, median, CV, "
                  f"p_traj = share of complete trajectories ({k1} points):")
            print(wide.round(3).to_string())
        print(f"\nReference: tip speed {0.576 * 10:.1f} mm/d. An output is usable for Morris if p_traj is close to 1.")


def compare_roots(out, design):
    """(on - off) / mean(off) per output and time, paired over the seeds, with 95% CI."""
    seed = design.set_index("run_id")[["group", "seed"]]
    d = out.join(seed, on="run_id")
    rows = []
    for (o, h), s in d[d.group.isin([OFF, ON])].groupby(["output", "time"]):
        w = s.pivot_table(index="seed", columns="group", values="value")
        if OFF not in w or ON not in w:
            continue
        diff = (w[ON] - w[OFF]).dropna()
        m = w[OFF].mean()
        if len(diff) >= 3 and m:
            rel = diff.mean() / m
            hw = stats.t.ppf(0.975, len(diff) - 1) * diff.std() / np.sqrt(len(diff)) / abs(m)
        else:
            rel = hw = np.nan
        rows.append({"output": o, "time": h, "n_pairs": len(diff), "mean_off": m, "mean_on": w[ON].mean(),
                     "rel_diff": rel, "rel_ci_lo": rel - hw, "rel_ci_hi": rel + hw})
    tab = pd.DataFrame(rows, columns=["output", "time", "n_pairs", "mean_off", "mean_on", "rel_diff", "rel_ci_lo", "rel_ci_hi"])
    tab.to_csv(config.RESULTS_DIR / "long_roots_compare.csv", index=False)
    if tab.empty:
        print("\nRoots on vs. off: no pairs (group long_roots_on not simulated yet?)")
        return
    print("\nRoots on vs. off: (on - off) / off, 95% CI paired over seeds  (<-- : CI outside ±5 %)")
    times = sorted(tab.time.unique())
    for o in [o for o in SHOW if o in set(tab.output)]:
        cells = []
        for h in times:
            r = tab[(tab.output == o) & (tab.time == h)]
            if r.empty or np.isnan(r.rel_diff.iloc[0]):
                cells.append(f"{h:g} h: -")
                continue
            r = r.iloc[0]
            flag = " <--" if (r.rel_ci_lo > 0.05 or r.rel_ci_hi < -0.05) else ""
            cells.append(f"{h:g} h: {r.rel_diff:+.3f} [{r.rel_ci_lo:+.3f}, {r.rel_ci_hi:+.3f}]{flag}")
        print(f"  {o:16s} " + " | ".join(cells))


def influx(design):
    """Hourly supply series of all runs with diagnostics, as a long table (hour after t*)."""
    rows, missing = [], []
    for _, r in design.iterrows():
        meta_path = config.RUNS_DIR / r.run_id / "meta.json"
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text())
        if meta.get("sim_status") != "ok":
            continue
        raw = model_hook.load_raw(config.RUNS_DIR / r.run_id / meta["snapshots"][0])
        if "diag_edges" not in raw:
            missing.append(r.run_id)
            continue
        hour = np.round((raw["diag_edges"][:-1] - float(raw["time_star"])) * 24).astype(int)
        df = pd.DataFrame({"run_id": r.run_id, "group": r.group, "seed": r.seed, "hour": hour})
        for k in INFLUX:
            df[k] = raw[k] if k in raw else np.nan
        rows.append(df)
    if missing:
        print(f"\n{len(missing)} runs without supply diagnostics (simulated before they existed), e.g. {missing[0]}. "
              "To include them: delete runs/<run_id>/ and resubmit the design.")
    if not rows:
        print("No runs with supply diagnostics.")
        return None
    per_run = pd.concat(rows, ignore_index=True)
    mean = per_run.groupby(["group", "hour"])[list(INFLUX)].mean().reset_index()
    n = per_run.groupby(["group", "hour"]).run_id.nunique().rename("n_runs").reset_index()
    mean = mean.merge(n, on=["group", "hour"])
    mean.to_csv(config.RESULTS_DIR / "long_influx.csv", index=False)

    print(f"\nSupply of hyphae to the fungal side (mean over seeds), crossings per hour:")
    for g, d in mean.groupby("group"):
        d = d.set_index("hour").sort_index()
        c = d["diag_cross_in"].rolling(6, min_periods=1).mean()
        after = c[c.index >= 0]
        windows = [(0, 10), (50, 60), (110, 120), (170, 180)]
        txt = ", ".join(f"{a}-{b} h: {d.loc[(d.index >= a) & (d.index < b), 'diag_cross_in'].mean():.1f}"
                        for a, b in windows if (d.index >= a).any() and (d.index < b).any())
        peak = c.max()
        low = after[after < DROP * peak]
        stop = f"below {DROP:.0%} of the maximum from {low.index[0]} h after t*" if len(low) else \
            f"not below {DROP:.0%} of the maximum until {after.index.max()} h after t*"
        print(f"  {LABELS.get(g, g)}: {txt}; 6 h mean {stop}")
    return mean


def plot_influx(mean):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    keys = [k for k in INFLUX if mean[k].notna().any()]
    fig, axes = plt.subplots(len(keys), 1, figsize=(8, 2.4 * len(keys) + 0.8), sharex=True, squeeze=False)
    for ax, k in zip(axes[:, 0], keys):
        for g, d in mean.groupby("group"):
            d = d.sort_values("hour")
            ax.plot(d.hour, d[k].rolling(3, center=True, min_periods=1).mean(), color=COLORS.get(g, INK2), lw=2)
        for x in (0, config.HOURS_HYPHAE):
            ax.axvline(x, color=GRID, lw=1, zorder=0)
        ax.set_title(INFLUX[k], fontsize=10, color=INK, loc="left")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
    axes[-1, 0].set_xlabel(f"h after t*  (grey lines: t* and {config.HOURS_HYPHAE} h, the window of SA and calibration)",
                           fontsize=8, color=INK2)
    handles = [Line2D([], [], color=COLORS[g], lw=2, label=LABELS[g]) for g in (OFF, ON) if g in set(mean.group)]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=9)
    fig.suptitle("Supply of hyphae, mean over seeds (3 h moving mean)", x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    (config.RESULTS_DIR / "plots").mkdir(parents=True, exist_ok=True)
    fig.savefig(config.RESULTS_DIR / "plots" / "long_influx.png", dpi=150)
    plt.close(fig)


def main():
    design = pd.read_csv(config.RESULTS_DIR / "design_long.csv")
    if "group" not in design:
        design["group"] = OFF
    status = pd.read_csv(config.RESULTS_DIR / "status_long.csv")
    out = pd.read_csv(config.RESULTS_DIR / "outputs_long.csv")
    summary(out, status, design)
    compare_roots(out, design)
    mean = influx(design)
    if mean is not None:
        plot_influx(mean)
    print(f"-> {config.RESULTS_DIR / 'long_summary.csv'}, long_roots_compare.csv, long_influx.csv, plots/long_influx.png")


if __name__ == "__main__":
    main()
