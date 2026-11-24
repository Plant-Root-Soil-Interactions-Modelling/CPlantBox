"""
Evaluates the raw output of ONE run -> runs/<run_id>/outputs.csv

    python evaluate.py results/design_morris.csv 17
    python evaluate.py results/design_morris.csv 0-9

Separate from the simulation, so that new outputs or times can be computed later without
re-simulating (then with --force).

An output is identified by (time, region, output); time in h after t*, region "all"
is the whole hyphal side of the dish. Two kinds:
  - quantities at every time in config.SNAPSHOT_HOURS (config.SNAPSHOT_OUTPUTS)
  - wave quantities over the whole window, at the last time (config.WAVE_OUTPUTS), from the
    fits of AMFDensities as in petridishstats.py
Plus diagnostics starting with "n_" (number of rings behind the fits), which do not enter the SA.
"""
import json
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

import config
import model_hook
from run_one import parse_rows

model_hook.import_amf()
import AMFDensities as amfd  # noqa: E402

FWHM_BELL = amfd.FWHM_BELL


def ring_centers_mm(location):
    inner = np.concatenate(([0.], location[:-1]))
    return np.sqrt((inner * inner + location * location) / 2) * 10   # amfs.ringCenters


def snapshot_outputs(raw, hour):
    j = int(np.argmin(np.abs(raw["times"] * 24 - hour)))      # column closest to the time
    area = raw["area"]
    L = {k: float((raw["length_" + k][:, j] * area).sum() / 10) for k in ("runner", "bas")}   # [cm]
    tips = raw["active_tips"][:, j]
    r = ring_centers_mm(raw["location"])
    return {
        "total_length": float((raw["length_densities"][:, j] * area).sum() / 10),   # [cm]
        "active_tips": float(tips.sum()),
        "bas_share": L["bas"] / (L["bas"] + L["runner"]) if L["bas"] + L["runner"] > 0 else np.nan,
        "anastomoses": float(raw["ana_new"][:, :j + 1].sum()),                       # since t*
        "front_position": float((r * tips).sum() / tips.sum()) if tips.sum() > 0 else np.nan,   # tip-weighted radius [mm]
        # without a fit, defined in every run:
        "spread_r90": spread_radius(raw, j, 0.9),                                    # 90 % of the length inside [mm]
        "max_tip_density": float(np.max(raw["tip_dens"][:, j])),                     # densest ring [mm^-2]
    }


def spread_radius(raw, j, q):
    # radius within which the share q of the hyphal length lies (linear within a ring) [mm]
    L = raw["length_densities"][:, j] * raw["area"].ravel()
    if L.sum() <= 0:
        return np.nan
    edges = np.concatenate(([0.], raw["location"])) * 10
    cum = np.concatenate(([0.], np.cumsum(L))) / L.sum()
    return float(np.interp(q, cum, edges))


def front_speed(raw, hour, start=10):
    # without a fit: slope of the tip-weighted radius over time (Theil-Sen) between start and hour [mm/d].
    # Can be negative (front retreats), NaN only without growing tips
    r = ring_centers_mm(raw["location"])
    t = raw["times"] * 24
    cols = np.where((t >= start - 1e-6) & (t <= hour + 1e-6))[0]
    tips = raw["active_tips"][:, cols]
    n = tips.sum(0)
    ok = n > 0
    if ok.sum() < 3:
        return np.nan
    fp = (r[:, None] * tips[:, ok]).sum(0) / n[ok]
    return float(amfd.theilSlope(t[cols][ok], fp)[0] * 24)


def truncate(raw, hour):
    # window up to hour [h after t*]
    keep = np.concatenate(([True], raw["times"][1:] * 24 <= hour + 1e-6))
    sim = {k: (v[:, keep] if isinstance(v, np.ndarray) and v.ndim == 2 else v) for k, v in raw.items()}
    sim["times"] = raw["times"][keep]
    return sim


def fit_table(sim):
    # wave fits per ring (AMFDensities.collectArrivalFits) as a table, sim from truncate
    r = ring_centers_mm(sim["location"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fits = amfd.collectArrivalFits([sim])
    rows = []
    for f in fits:
        K1, lam_len, t_n = f["length"]
        K2, lam_tip, t_tip = f["tips"] if f["tips"] is not None else (np.nan, np.nan, np.nan)
        rows.append({"r": r[f["ring"]], "t_n": t_n, "K1": K1, "R2_length": f["R2_length"],
                     "saturation": f["rho_data"][-1] / K1, "K2": K2, "lam_tip": lam_tip, "t_tip": t_tip,
                     "R2_tips": f["R2_tips"], "t_first": f["t"][0], "t_last": f["t"][-1]})
    return pd.DataFrame(rows, columns=["r", "t_n", "K1", "R2_length", "saturation", "K2", "lam_tip", "t_tip",
                                       "R2_tips", "t_first", "t_last"])


def good_rings(df):
    # quality criteria as petridishstats.py -> (all good, saturated ones for K1, tip peak in the window for K2)
    good = df[df.R2_length >= config.R2_MIN_LENGTH]
    good_K1 = good[good.saturation >= config.SATURATION_MIN]
    # as petridishwave.py: tip peak and both half-maximum points inside the simulated window, otherwise extrapolated
    half = FWHM_BELL / good.lam_tip / 2
    inside = (good.t_tip - half >= good.t_first) & (good.t_tip + half <= good.t_last)
    good_tips = good[(good.R2_tips >= config.R2_MIN_TIPS) & inside]
    return good, good_K1, good_tips


def wave_outputs(raw, hour):
    sim = truncate(raw, hour)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        C = amfd.encounterC([sim], min_tips=config.ENCOUNTER_MIN_TIPS)["total"]
        beta = amfd.anastomosisBeta([sim], min_tips=config.ENCOUNTER_MIN_TIPS)["total"]
    good, good_K1, good_tips = good_rings(fit_table(sim))

    def velocity(t, rr):
        # as petridishwave.py: radius vs. arrival time (Theil-Sen) [mm/h] -> [mm/d]; rising fronts only
        if len(t) < config.MIN_RINGS:
            return np.nan
        s = amfd.theilSlope(t, rr)[0]
        return 24 * s if np.isfinite(s) and s > 0 else np.nan

    med = lambda d, x: float(np.median(x)) if len(d) >= config.MIN_RINGS else np.nan
    return {
        "front_speed": front_speed(raw, hour),                          # [mm/d] without a fit, from front_position(t)
        "front_velocity": velocity(good_tips.t_tip, good_tips.r),      # [mm/d] from the tip peak per ring (v_tips)
        # comparison variants (not SA outputs): arrival t_n of the length density in all good or saturated rings
        "v_length": velocity(good.t_n, good.r),
        "v_length_sat": velocity(good_K1.t_n, good_K1.r),
        "K1": med(good_K1, good_K1.K1),                                 # [µm mm^-2]
        "K2": med(good_tips, good_tips.K2),                             # [mm^-2]
        "wave_width": med(good_tips, FWHM_BELL / good_tips.lam_tip),    # duration of the tip wave in a ring [h]
        "encounter_C": float(C),                                        # f = C rho [mm^2 µm^-1 h^-1]
        "beta": float(beta),                                            # a = beta rho n [µm h^-1], target ≈ 23 (Oyarte Galvez 2025)
        "n_wave_rings": len(good), "n_K1_rings": len(good_K1), "n_tip_rings": len(good_tips),
    }


def evaluate_run(run_id, force=False):
    rdir = config.RUNS_DIR / run_id
    meta_path = rdir / "meta.json"
    if not meta_path.exists():
        print(f"{run_id}: not simulated yet"); return
    meta = json.loads(meta_path.read_text())
    out_path = rdir / "outputs.csv"
    if meta.get("sim_status") != "ok":
        print(f"{run_id}: simulation failed, no evaluation"); return
    if out_path.exists() and not force:
        print(f"{run_id}: already evaluated"); return
    raw = model_hook.load_raw(rdir / meta["snapshots"][0])
    rows = []
    simulated = raw["times"][-1] * 24 + 1                              # [h after t*], +1 h tolerance
    extra = [h for h in config.EXTRA_HOURS if max(config.SNAPSHOT_HOURS) < h <= simulated]   # long runs only
    for h in config.SNAPSHOT_HOURS + extra:
        for k, v in snapshot_outputs(raw, h).items():
            rows.append({"run_id": run_id, "time": h, "region": "all", "output": k, "value": v})
    for h in [max(config.SNAPSHOT_HOURS)] + extra:                     # wave fits over the window 0..h
        for k, v in wave_outputs(raw, h).items():
            rows.append({"run_id": run_id, "time": h, "region": "all", "output": k, "value": v})
    rows.append({"run_id": run_id, "time": 0, "region": "all", "output": "wait_days", "value": float(raw["wait_days"])})
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f"{run_id}: {len(rows)} values")


if __name__ == "__main__":
    design = pd.read_csv(sys.argv[1])
    force = "--force" in sys.argv
    for i in parse_rows(sys.argv[2]):
        if i < len(design):
            evaluate_run(design.iloc[i]["run_id"], force)
