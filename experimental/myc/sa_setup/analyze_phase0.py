"""
Phase 0: seed noise and dx convergence at the nominal point.

    python analyze_phase0.py
      -> results/phase0_noise.csv   mean, SD, CV, share of zeros per output
                                    + noise level of the elementary effects for R = 1..5
      -> results/phase0_dx.csv      relative difference per dx vs. the finest dx
      -> results/phase0_dt.csv      relative difference per time step dt vs. the finest dt
      -> results/phase0_status.csv  share of runs with t* and runtime per group (range: sample of the factor space)
"""
import numpy as np
import pandas as pd
from scipy import stats

import config

KEY = ["time", "region", "output"]


def convergence(out, var):
    """Relative difference of every output between each value of var (groups "dx=0.005", "dt=0.5h") and the finest,
    paired over the seeds, with 95% CI. Empty table if the design has no such groups."""
    cols = list(KEY) + [var, "ref_" + var, "n_pairs", "rel_diff", "rel_ci_lo", "rel_ci_hi"]
    sub_all = out[out.group.str.startswith(var + "=")].copy()
    if sub_all.empty:
        return pd.DataFrame(columns=cols)
    sub_all[var] = sub_all.group.str[len(var) + 1:].str.rstrip("h").astype(float)
    finest = sub_all[var].min()
    rows = []
    for key, sub in sub_all.groupby(KEY):
        ref = sub[sub[var] == finest].set_index("seed")["value"]
        for val, s in sub.groupby(var):
            if val == finest:
                continue
            diff = (s.set_index("seed")["value"] - ref).dropna()
            m = ref.mean()
            if len(diff) >= 3 and m:      # fewer pairs (e.g. output often NaN): no reliable comparison
                rel = diff.mean() / m
                hw = stats.t.ppf(0.975, len(diff) - 1) * diff.std() / np.sqrt(len(diff)) / abs(m)
            else:
                rel = hw = np.nan
            rows.append(dict(zip(KEY, key), **{var: val, "ref_" + var: finest}, n_pairs=len(diff), rel_diff=rel,
                             rel_ci_lo=rel - hw, rel_ci_hi=rel + hw))
    return pd.DataFrame(rows, columns=cols)


def print_convergence(tab, var, nominal, unit, t_end):
    d = tab[(tab.time == t_end) & (tab.region == "all") & (tab[var] == nominal)]
    if d.empty:
        print(f"\n{var} convergence: no runs")
        return
    print(f"\n{var} convergence: {var} = {nominal:g}{unit} vs. {var} = {d['ref_' + var].iloc[0]:g}{unit} "
          "(rel. difference, 95% CI):")
    for _, r in d.iterrows():
        if np.isnan(r.rel_diff):
            print(f"  {r.output:16s} too few pairs ({r.n_pairs})")
            continue
        flag = "  <-- check" if (r.rel_ci_lo > 0.02 or r.rel_ci_hi < -0.02) else ""
        print(f"  {r.output:16s} {r.rel_diff:+.3f} [{r.rel_ci_lo:+.3f}, {r.rel_ci_hi:+.3f}]  n={r.n_pairs}{flag}")


def main():
    design = pd.read_csv(config.RESULTS_DIR / "design_phase0.csv")
    out = pd.read_csv(config.RESULTS_DIR / "outputs_phase0.csv").merge(
        design[["run_id", "group", "seed"]], on="run_id")

    # ---------- noise ----------
    noise = out[out.group == "noise"]
    g = noise.groupby(KEY)["value"]
    tab = pd.DataFrame({"n": g.size(), "n_valid": g.count(), "mean": g.mean(), "median": g.median(), "sd": g.std(),
                        "frac_zero": g.apply(lambda v: np.mean(v == 0))})
    tab["cv"] = tab["sd"] / tab["mean"].replace(0, np.nan)
    # noise level of an elementary effect (normalised units), averaged over R replicates,
    # conservatively without the gain from shared seeds: sqrt(2)*sd / (Δ*sqrt(R)).
    p = config.MORRIS_LEVELS
    delta = p / (2 * (p - 1))
    for R in range(1, 6):
        tab[f"EE_noise_R{R}"] = np.sqrt(2) * tab["sd"] / (delta * np.sqrt(R))
    tab.reset_index().to_csv(config.RESULTS_DIR / "phase0_noise.csv", index=False)

    # ---------- dx and dt convergence (paired over seeds) ----------
    dxt = convergence(out, "dx")
    dxt.to_csv(config.RESULTS_DIR / "phase0_dx.csv", index=False)
    dtt = convergence(out, "dt")
    dtt.to_csv(config.RESULTS_DIR / "phase0_dt.csv", index=False)

    # ---------- sample of the factor space: aborts and runtime ----------
    status = pd.read_csv(config.RESULTS_DIR / "status_phase0.csv").merge(design[["run_id", "group"]], on="run_id")
    for g in ("dx", "dt"):
        status["group"] = status["group"].where(~status["group"].str.startswith(g + "="), g)
    rng_tab = status.groupby("group").agg(n=("status", "size"), ok=("status", lambda s: int((s == "ok").sum())),
                                          sim_h_median=("sim_seconds", lambda s: s.median() / 3600),
                                          sim_h_max=("sim_seconds", lambda s: s.max() / 3600))
    rng_tab.to_csv(config.RESULTS_DIR / "phase0_status.csv")

    # ---------- short report ----------
    t_end = max(config.SNAPSHOT_HOURS)
    print("\nRuns per group (ok = t* reached and evaluated), simulation runtime [h]:")
    print(rng_tab.round(2).to_string())
    print(f"\nSeed noise at the nominal point (t = {t_end} h after t*):")
    print(tab.loc[(t_end, "all")][["n_valid", "mean", "median", "sd", "cv", "frac_zero"]].round(3).to_string())
    print_convergence(dxt, "dx", max(config.PHASE0_DX_VALUES), "", t_end)
    print_convergence(dtt, "dt", 24 / config.STEPS_PER_DAY, " h", t_end)
    sparse = tab[tab.frac_zero > 0.2].reset_index()
    if len(sparse):
        print(f"\n{len(sparse)} output/region/time combinations with >20 % zeros "
              "(sparsely populated subregions, see phase0_noise.csv)")


if __name__ == "__main__":
    main()
