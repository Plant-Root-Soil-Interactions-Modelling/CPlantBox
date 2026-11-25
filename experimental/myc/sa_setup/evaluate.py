"""
Evaluates the raw output of ONE run -> runs/<run_id>/outputs.csv

    python evaluate.py results/design_morris.csv 17
    python evaluate.py results/design_morris.csv 0-9

Separate from the simulation, so that new outputs or times can be computed later without
re-simulating (then with --force).

An output is identified by (time, region, output); time in h after t*, region "all"
is the whole hyphal side of the dish. Two kinds:
  - quantities at every time in config.SNAPSHOT_HOURS (config.SNAPSHOT_OUTPUTS)
  - wave quantities over the whole window, at the last time (WAVE_QUANTITIES, one definition for SA and calibration),
    from the fits of AMFDensities as in petridishstats.py; pooled per design point: evaluate_pooled.py
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
    n_ana = float(raw["ana_new"][:, :j + 1].sum())                                   # anastomoses since t*
    return {
        "total_length": float((raw["length_densities"][:, j] * area).sum() / 10),   # [cm]
        "active_tips": float(tips.sum()),
        "bas_share": L["bas"] / (L["bas"] + L["runner"]) if L["bas"] + L["runner"] > 0 else np.nan,
        "anastomoses": n_ana,
        "front_position": float((r * tips).sum() / tips.sum()) if tips.sum() > 0 else np.nan,   # tip-weighted radius [mm]
        # without a fit, defined in every run:
        "spread_r90": spread_radius(raw, j, 0.9),                                    # 90 % of the length inside [mm]
        "max_tip_density": float(np.max(raw["tip_dens"][:, j])),                     # densest ring [mm^-2]
        # image graph (only at config.IMAGE_HOURS)
        **image_outputs(raw, hour, n_ana),
        # diagnostics (E3, no SA outputs, no targets): where do the anastomoses happen, how many growing tips are BAS?
        **bas_diagnostics(raw, j),
    }


def bas_diagnostics(raw, j):
    """Share of the anastomoses since t* made by BAS tips, share of BAS among the growing tips at column j (E3)."""
    out = {}
    if "ana_new_bas" in raw:
        a_all = float(raw["ana_new"][:, :j + 1].sum())
        out["ana_share_bas"] = float(raw["ana_new_bas"][:, :j + 1].sum()) / a_all if a_all > 0 else np.nan
    if "active_tips_bas" in raw:
        t_all = float(raw["active_tips"][:, j].sum())
        out["tip_share_bas"] = float(raw["active_tips_bas"][:, j].sum()) / t_all if t_all > 0 else np.nan
    return out


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


# --------------------------------------------------------------------------
# Wave quantities: ONE definition for SA, calibration and all analyses (B12)
# --------------------------------------------------------------------------
# Nomenclature: alpha = anastomosis (a = alpha rho n), beta = branching (b = beta n). Oyarte Galvez et al. (2025) use the
# letters the other way round. As in the experiment (Methods: only tips that moved count, BAS tips hardly), tips are the
# growing RUNNER tips: tip density of the wave fits, n in alpha, beta. rho is the filament density of all hyphae.
WAVE_QUANTITIES = {
    "alpha": "[µm/h] anastomosis a = alpha rho n, a and n of the runner tips (ratio of sums), target ≈ 23",
    "beta": "[h^-1] branching: new runner tips per growing runner tip and hour, target ≈ 0.04",
    "K1": "[µm mm^-2] saturation of the filament density (median of saturated rings)",
    "K2": "[mm^-2] peak of the runner tip density (median of rings with the peak in the window)",
    "wave_width": "[h] duration of the runner tip wave in a ring (median)",
    "front_velocity": "[mm/d] ring radius vs. arrival t_n of the filament density (sigmoid fit, Theil-Sen)",
    "front_speed": "[mm/d] slope of the tip-weighted radius from 10 h, without a fit (mean over runs)",
    "encounter_C": "[mm^2 µm^-1 h^-1] f = C rho, all tips",
    "front_tip_speed": "[µm/h] growth speed of the runner tips at the front (outermost FRONT_QUANTILE), mean over the hours "
                       "from SPEED_FROM (Fig. 2d: ≈ 300)",
    "tip_speed_mean": "[µm/h] growth speed of all growing runner tips on the fungal side, same hours (comparison: v ≈ 240)",
    "branch_angle_rh": "[°] mean absolute branching angle of new-born runners, x-y projection, initial growth direction "
                       "(Extended Data Fig. 11b: 73°), end of the simulation",
    "branch_angle_bas": "[°] same for new-born BAS (from runners and tip splits) (80°)",
    "persistence_length": "[cm] median persistence length of the runners, <cos Δθ(ΔL)> = exp(-ΔL/l_p) up to 1 cm "
                          "(Extended Data Fig. 11a: ≈ 10 cm, a lower bound)",
}


def _speed_mean(sim, key, hour):
    """Mean of the hourly tip speed `key` (model_hook.tip_speed_diagnostics) over SPEED_FROM < h <= hour, nan without it."""
    if key not in sim:
        return np.nan
    h, v = np.asarray(sim["diag_speed_hours"]), np.asarray(sim[key], dtype=float)
    sel = (h > config.SPEED_FROM - 1e-6) & (h <= hour + 1e-6)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(np.nanmean(v[sel])) if sel.any() else np.nan

# Quantities of the planar image graph at config.IMAGE_HOURS after t* (model_hook, from the final plant), per run like the
# snapshot outputs (calibration: mean over the seeds at Target.hour)
IMAGE_QUANTITIES = {
    "root_efficiency_rel": "(E_r - E_r,MST) / (E_r,DT - E_r,MST), 2D (figure: Ê_r ≈ 0.69, 0.665, 0.63 at 20, 40, 60 h)",
    "bas_share_image": "BAS share of the length as the experiment classifies the edges of the image graph (edge length, "
                       "tip, width; AMFNetwork.basShareImage) (ρ_BAS ≈ 0.37, 0.385, 0.40)",
    "bas_share_subtype_image": "comparison: BAS share by sub type in the same image graph",
    "L_RH_image": "[mm] runner length as classified in the image graph, projected (≈ 280, 500, 850)",
    "L_BAS_image": "[mm] BAS length as classified in the image graph, projected (≈ 100, 230, 370)",
    # anastomosis density ρ_A of the paper (≈ 0.5–0.7, unit and normalisation to be checked): both candidates
    "ana_per_length_rh": "[mm^-1] anastomoses since t* (fungal side) per L_RH of the image graph",
    "ana_per_length_total": "[mm^-1] anastomoses since t* (fungal side) per total length of the image graph",
    # held back for the validation (config.VALIDATION_OUTPUTS), normalisations assumed
    "global_efficiency_rel": "(E_g - E_g,MST) / (E_g,DT - E_g,MST) over the pairs of terminal and root nodes (Ê_g ≈ 0.33, 0.44, 0.53)",
    "cost_rel": "(L - L_MST) / (L_DT - L_MST), total length of the image graph vs. the reference networks (Ĉ ≈ 0.1)",
}


def _image_value(raw, key, hour):
    """raw[<key>_<hour>h] (model_hook, IMAGE_HOURS); at the end of the run also the key without suffix (older raw.npz)."""
    k = f"{key}_{hour:g}h"
    if k in raw:
        return float(raw[k])
    if key in raw and abs(raw["times"][-1] * 24 - hour) < 1.5:
        return float(raw[key])
    return None


def image_outputs(raw, hour, n_ana):
    """IMAGE_QUANTITIES at `hour` after t*, {} if the image graph was not computed at that time. n_ana: anastomoses since t*."""
    share = _image_value(raw, "bas_image_share", hour)
    if share is None:
        return {}
    L = _image_value(raw, "bas_image_length_xy_um", hour)
    L = np.nan if L is None else L / 1000.                       # [mm]
    sub = _image_value(raw, "bas_image_share_subtype", hour)
    eff = _image_value(raw, "root_eff_E_rel", hour)
    eg = _image_value(raw, "root_eff_E_g_rel", hour)
    cost = _image_value(raw, "root_eff_C_rel", hour)
    L_rh, L_bas = (1 - share) * L, share * L
    with np.errstate(invalid="ignore", divide="ignore"):
        return {"bas_share_image": share, "bas_share_subtype_image": np.nan if sub is None else sub,
                "L_RH_image": L_rh, "L_BAS_image": L_bas,
                "root_efficiency_rel": np.nan if eff is None else eff,
                "global_efficiency_rel": np.nan if eg is None else eg, "cost_rel": np.nan if cost is None else cost,
                "ana_per_length_rh": n_ana / L_rh if L_rh > 0 else np.nan,
                "ana_per_length_total": n_ana / L if L > 0 else np.nan}


def prepare(sim, hour):
    """Expensive part per run (fits), sim from truncate(raw, hour). The pooled quantities reuse it (bootstrap)."""
    sim_rh = dict(sim, tip_dens=sim["active_tips_runner"] / np.asarray(sim["area"]).reshape(-1, 1))
    return {"sim": sim, "table": fit_table(sim_rh), "front_speed": front_speed(sim, hour),
            "front_tip_speed": _speed_mean(sim, "diag_front_tip_speed", hour),
            "tip_speed_mean": _speed_mean(sim, "diag_tip_speed_all", hour)}


def wave_quantities(preps):
    """WAVE_QUANTITIES pooled over the runs of `preps` (one run: per-run value), plus comparison variants and the number
    of rings behind the fits (n_*). Used per run (evaluate), per design point (evaluate_pooled) and in the calibration."""
    sims = [p["sim"] for p in preps]
    good, good_K1, good_tips = good_rings(pd.concat([p["table"] for p in preps], ignore_index=True))
    n_min = config.MIN_RINGS
    med = lambda d, x: float(np.median(x)) if len(d) >= n_min else np.nan

    def velocity(t, rr):
        return 24 * float(amfd.theilSlope(t, rr)[0]) if len(t) >= n_min else np.nan

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return {
            "alpha": float(amfd.anastomosisAlpha(sims, min_tips=config.ENCOUNTER_MIN_TIPS)["runner"]),
            "beta": float(amfd.branchingRate(sims, min_tips=config.ENCOUNTER_MIN_TIPS, suffix="_runner")),
            "K1": med(good_K1, good_K1.K1),
            "K2": med(good_tips, good_tips.K2),
            "wave_width": med(good_tips, FWHM_BELL / good_tips.lam_tip),
            "front_velocity": velocity(good.t_n, good.r),
            "front_speed": float(np.nanmean([p["front_speed"] for p in preps])),
            "front_tip_speed": float(np.nanmean([p["front_tip_speed"] for p in preps])),
            "tip_speed_mean": float(np.nanmean([p["tip_speed_mean"] for p in preps])),
            "branch_angle_rh": float(np.nanmean([float(s.get("angle_rh_mean", np.nan)) for s in sims])),
            "branch_angle_bas": float(np.nanmean([float(s.get("angle_bas_mean", np.nan)) for s in sims])),
            "persistence_length": float(np.nanmedian([float(s.get("persistence_length_median", np.nan)) for s in sims])),
            "encounter_C": float(amfd.encounterC(sims, min_tips=config.ENCOUNTER_MIN_TIPS)["total"]),
            # comparison variants (not SA outputs or targets): from the runner tip peak, from saturated rings only
            "v_tips": velocity(good_tips.t_tip, good_tips.r),
            "v_length_sat": velocity(good_K1.t_n, good_K1.r),
            "n_wave_rings": len(good), "n_K1_rings": len(good_K1), "n_tip_rings": len(good_tips),
        }


def wave_outputs(raw, hour):
    """Wave quantities of ONE run over the window 0..hour after t*."""
    sim = truncate(raw, hour)
    return wave_quantities([prepare(sim, hour)])


def influx_outputs(raw, hour):
    """Supply of hyphae to the fungal side from model_hook.influx_diagnostics (crossings per hour), {} without them:
    influx_tstar: mean over the 6 h before t*, influx_window: mean over 0..hour after t*,
    influx_decay: mean over the last 10 h of the window / influx_tstar (1 = constant supply, 0 = supply stopped)."""
    if "diag_edges" not in raw:
        return {}
    h = (raw["diag_edges"][:-1] - float(raw["time_star"])) * 24
    c = raw["diag_cross_in"]
    mean = lambda a, b: float(c[(h >= a - 1e-6) & (h < b - 1e-6)].mean()) if ((h >= a - 1e-6) & (h < b - 1e-6)).any() else np.nan
    at_star = mean(-6, 0)
    out = {"influx_tstar": at_star, "influx_window": mean(0, hour),
           "influx_decay": mean(hour - 10, hour) / at_star if at_star > 0 else np.nan}
    if "diag_frozen_tips_fungal" in raw:
        # share of the tips on the fungal side that are frozen because an ancestor died (artefact of Hyphae::simulate, A6),
        # at the end of the window: frozen / (frozen + growing)
        k = int(np.argmin(np.abs(h + 1 - hour)))      # bin ending at `hour`
        frozen = float(raw["diag_frozen_tips_fungal"][k])
        j = int(np.argmin(np.abs(raw["times"] * 24 - hour)))
        growing = float(raw["active_tips"][:, j].sum())
        out["frozen_share"] = frozen / (frozen + growing) if frozen + growing > 0 else np.nan
        out["frozen_tips"] = frozen
    return out


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
        for k, v in {**wave_outputs(raw, h), **influx_outputs(raw, h)}.items():
            rows.append({"run_id": run_id, "time": h, "region": "all", "output": k, "value": v})
    rows.append({"run_id": run_id, "time": 0, "region": "all", "output": "wait_days", "value": float(raw["wait_days"])})
    for k in ("ring_center_y", "n_crossings_tstar"):   # barrier crossing point of the run (centre of the rings), diagnostics
        if k in raw:
            rows.append({"run_id": run_id, "time": 0, "region": "all", "output": k, "value": float(raw[k])})
    if "diag_frac_tips_outside_rings" in raw:          # share of the growing runner tips beyond RING_RADIUS, diagnostics
        h, f = np.asarray(raw["diag_speed_hours"]), np.asarray(raw["diag_frac_tips_outside_rings"], dtype=float)
        for t in config.SNAPSHOT_HOURS + extra:
            j = np.argmin(np.abs(h - t))
            rows.append({"run_id": run_id, "time": t, "region": "all", "output": "frac_tips_outside_rings", "value": float(f[j])})
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f"{run_id}: {len(rows)} values")


if __name__ == "__main__":
    design = pd.read_csv(sys.argv[1])
    force = "--force" in sys.argv
    for i in parse_rows(sys.argv[2]):
        if i < len(design):
            evaluate_run(design.iloc[i]["run_id"], force)
