"""
Phase 0: seed noise and dx convergence at the nominal point.

    python analyze_phase0.py
      -> results/phase0_noise.csv   mean, SD, CV, share of zeros per output
                                    + noise level of the elementary effects for R = 1..5
      -> results/phase0_dx.csv      relative difference per dx vs. the finest dx
      -> results/phase0_dt.csv      relative difference per time step dt vs. the finest dt
      -> results/phase0_layer.csv   relative difference per thickness of the hyphal layer vs. the thinnest
      -> results/phase0_anaL.csv    relative difference per factor on the derived anaMinDist vs. the nominal factor 3
      -> results/phase0_status.csv  share of runs with t* and runtime per group (range: sample of the factor space)
"""
import numpy as np
import pandas as pd
from scipy import stats

import config

KEY = ["time", "region", "output"]


def convergence(out, var, ref=None):
    """Relative difference of every output between each value of var (groups "dx=0.005", "dt=0.5h") and the reference
    (default: the finest = smallest value), paired over the seeds, with 95% CI. Empty table if the design has no such groups."""
    cols = list(KEY) + [var, "ref_" + var, "n_pairs", "rel_diff", "rel_ci_lo", "rel_ci_hi"]
    sub_all = out[out.group.str.startswith(var + "=")].copy()
    if sub_all.empty:
        return pd.DataFrame(columns=cols)
    sub_all[var] = sub_all.group.str[len(var) + 1:].str.rstrip("h").astype(float)
    finest = sub_all[var].min() if ref is None else ref
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


def single_colony_cv(tab, t_end):
    """CV between the seeds at the nominal point of the targets that come from a single colony (calibrate.Target with
    replicates=1): the estimate of the spread between plates, to be entered by hand as Target.between_cv."""
    try:
        import calibrate
    except Exception:
        return
    single = [t for t in calibrate.TARGETS if t.replicates == 1]
    if not single:
        return
    print(f"\nSpread between seeds at the nominal point (estimate of the spread between plates) for the targets from a single"
          f"\ncolony -> enter as between_cv in calibrate.TARGETS (per run; the calibration pools, these are per-run CVs):")
    for t in single:
        key = (t.hour if t.hour is not None else t_end, "all", t.output)
        cv = tab["cv"].get(key, np.nan) if hasattr(tab["cv"], "get") else np.nan
        print(f"  {t.label:18s} CV {cv:.3f}   (now between_cv = {t.between_cv:g})")


def crossing_point(noise):
    """Barrier crossing point of the runs (centre of the rings, model_hook): distance from the centre of the dish."""
    s = noise[noise.output == "ring_center_y"].value
    n = noise[noise.output == "n_crossings_tstar"].value
    if len(s):
        print(f"\nBarrier crossing point (centre of the rings), y [cm] over the seeds: median {s.median():+.2f}, "
              f"range [{s.min():+.2f}, {s.max():+.2f}], crossings until t*: median {n.median():.0f}")


def ring_warnings(n_cross, frac_outside, where=""):
    """Lines "ATTENTION ..." for the rings (also used by calibrate.write_report): few barrier crossings until t* (uncertain
    crossing point = centre of the rings) and growing tips beyond RING_RADIUS (the rings miss the front)."""
    L = []
    n_cross, frac_outside = pd.Series(n_cross, dtype=float).dropna(), pd.Series(frac_outside, dtype=float).dropna()
    if len(n_cross) and n_cross.median() < config.WARN_MIN_CROSSINGS:
        L.append(f"ATTENTION number of barrier crossings until t*{where}: median {n_cross.median():.0f} "
                 f"(< {config.WARN_MIN_CROSSINGS}, min {n_cross.min():.0f}) -> the crossing point (centre of the rings) is "
                 "uncertain; average the crossings over a longer time (e.g. until 10 h after t*)")
    if len(n_cross) and (n_cross == 0).any():
        L.append(f"ATTENTION no barrier crossing until t* in {int((n_cross == 0).sum())} of {len(n_cross)} runs{where} "
                 "-> rings around the centre of the dish there")
    if len(frac_outside) and frac_outside.median() > config.WARN_MAX_OUTSIDE_RINGS:
        L.append(f"ATTENTION growing tips beyond RING_RADIUS = {config.RING_RADIUS * 10:g} mm{where}: median share "
                 f"{frac_outside.median():.0%} at the end (> {config.WARN_MAX_OUTSIDE_RINGS:.0%}, max {frac_outside.max():.0%}) "
                 "-> the front leaves the rings, K1 / K2 / wave fits miss the outer part; larger RING_RADIUS?")
    return L


def write_warnings(noise, t_end):
    n_cross = noise[noise.output == "n_crossings_tstar"].value
    outside = noise[(noise.output == "frac_tips_outside_rings") & (noise.time == t_end)].value
    L = ring_warnings(n_cross, outside, " (Phase 0 noise runs, nominal point)")
    if n_cross.empty or outside.empty:
        # runs without the ring diagnostics (simulated with an older version): no statement possible
        L.insert(0, "ATTENTION no ring diagnostics in the runs (n_crossings_tstar / frac_tips_outside_rings missing) -> "
                    "simulated with an older version of model_hook? The checks for the rings could not be made")
    path = config.RESULTS_DIR / "phase0_warnings.txt"
    path.write_text("\n".join(L) + "\n" if L else "no warnings\n", encoding="utf-8")
    print("\n" + ("\n".join(L) if L else "No warnings for the rings.") + f"\n-> {path}")


def crn(out, t_end):
    """Common random numbers (C19): sd of the difference moved − nominal with the same seeds (crn_same vs. the noise run of
    the same seed) and with independent seeds (crn_indep vs. the noise run with the same index). ratio < 1: shared seeds
    remove noise (ratio² = variance left); ~1: they bring nothing, the trajectories decorrelate."""
    nom = out[out.group == "noise"]
    same = out[out.group == "crn_same"]
    indep = out[out.group == "crn_indep"]
    if same.empty or indep.empty:
        return
    rows = []
    for (h, o), s in same.groupby(["time", "output"]):
        ref = nom[(nom.time == h) & (nom.output == o)].sort_values("seed")
        ref_by_seed = ref.set_index("seed").value
        d_same = (s.set_index("seed").value - ref_by_seed).dropna()
        ind = indep[(indep.time == h) & (indep.output == o)].sort_values("seed").value.values
        d_ind = pd.Series(ind[:len(ref)] - ref.value.values[:len(ind)]).dropna()
        if len(d_same) >= 3 and len(d_ind) >= 3 and d_ind.std() > 0:
            rows.append({"time": h, "output": o, "n": len(d_same), "mean_diff": d_same.mean(), "sd_same": d_same.std(),
                         "sd_indep": d_ind.std(), "ratio": d_same.std() / d_ind.std()})
    tab = pd.DataFrame(rows)
    tab.to_csv(config.RESULTS_DIR / "phase0_crn.csv", index=False)
    sub = tab[tab.time == t_end]
    if len(sub):
        print(f"\nCommon random numbers ({config.PHASE0_CRN_PARAM} moved by {config.PHASE0_CRN_STEP} normalised, {t_end} h):"
              "\n  sd(difference, same seeds) / sd(difference, independent seeds); < 1: shared seeds help")
        for _, r in sub.iterrows():
            print(f"  {r.output:22s} ratio {r.ratio:.2f}  (mean difference {r.mean_diff:+.4g}, n = {r.n})")


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
    # layer thickness: reference is the thinnest layer (closest to the 2D cellophane of the experiment)
    lyt = convergence(out, "layer")
    lyt.to_csv(config.RESULTS_DIR / "phase0_layer.csv", index=False)
    # minimal path length for anastomosis: factor on the derived value, reference = the nominal scale
    alt = convergence(out, "anaL", ref=config.ANA_L_SCALE_NOMINAL)
    alt.to_csv(config.RESULTS_DIR / "phase0_anaL.csv", index=False)

    # ---------- sample of the factor space: aborts and runtime ----------
    status = pd.read_csv(config.RESULTS_DIR / "status_phase0.csv").merge(design[["run_id", "group"]], on="run_id")
    for g in ("dx", "dt", "layer", "anaL"):
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
    print_convergence(lyt, "layer", config.HYPHAE_LAYER, " cm", t_end)
    for s in [s for s in config.PHASE0_ANAL_SCALES if s != config.ANA_L_SCALE_NOMINAL]:
        print_convergence(alt, "anaL", s, " × derived", t_end)
    # BAS share: by sub type (3D length, the snapshot output), by sub type in the image graph, and as the experiment
    # classifies the edges of the image graph (calibration target, Oyarte Galvez et al. 2025: 0.3)
    keys = [("bas_share", "by sub type, 3D length"), ("bas_share_subtype_image", "by sub type, image graph (2D)"),
            ("bas_share_image", "image criteria as the experiment (≈ 0.40 at 60 h)")]
    rows = []
    t_bas = max(config.IMAGE_HOURS)   # time of the experimental value
    for k, label in keys:
        v = noise[(noise.output == k) & (noise.time == t_bas)].set_index("seed").value
        rows.append((label, v))
    if any(len(v.dropna()) for _, v in rows):
        print(f"\nBAS share at the nominal point ({t_bas} h after t*), mean ± sd over the seeds:")
        for label, v in rows:
            print(f"  {label:48s} {v.mean():.3f} ± {v.std():.3f}  (n = {v.notna().sum()})")
        a, b = rows[0][1], rows[2][1]
        d = (b - a).dropna()
        if len(d) >= 3:
            print(f"  image criteria − sub type (paired): {d.mean():+.3f} ± {d.std():.3f}")
        cmp = pd.DataFrame({label: v for label, v in rows})
        cmp.to_csv(config.RESULTS_DIR / "phase0_bas_share.csv")
    crn(out, t_end)
    single_colony_cv(tab, t_end)
    crossing_point(noise)
    write_warnings(noise, t_end)
    sparse = tab[tab.frac_zero > 0.2].reset_index()
    if len(sparse):
        print(f"\n{len(sparse)} output/region/time combinations with >20 % zeros "
              "(sparsely populated subregions, see phase0_noise.csv)")


if __name__ == "__main__":
    main()
