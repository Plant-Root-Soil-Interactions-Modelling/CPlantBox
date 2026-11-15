"""
Wertet die Rohausgabe EINES Laufs aus -> runs/<run_id>/outputs.csv

    python evaluate.py results/design_morris.csv 17
    python evaluate.py results/design_morris.csv 0-9

Getrennt von der Simulation, damit neue Outputs oder Zeitpunkte später ohne
Neusimulation berechnet werden können (dann mit --force).

Ein Output wird über (time, region, output) identifiziert; time in h nach t*, region "all"
ist die ganze Hyphenseite der Schale. Zwei Arten:
  - Mengen zu jedem Zeitpunkt in config.SNAPSHOT_HOURS (config.SNAPSHOT_OUTPUTS)
  - Wellengrößen über das ganze Fenster, beim letzten Zeitpunkt (config.WAVE_OUTPUTS), aus den
    Fits von AMFDensities wie in petridishstats.py
Dazu Diagnosen mit "n_" (Anzahl Ringe hinter den Fits), die nicht in die SA eingehen.
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
    j = int(np.argmin(np.abs(raw["times"] * 24 - hour)))      # Spalte am nächsten zum Zeitpunkt
    area = raw["area"]
    L = {k: float((raw["length_" + k][:, j] * area).sum() / 10) for k in ("runner", "bas")}   # [cm]
    tips = raw["active_tips"][:, j]
    r = ring_centers_mm(raw["location"])
    return {
        "total_length": float((raw["length_densities"][:, j] * area).sum() / 10),   # [cm]
        "active_tips": float(tips.sum()),
        "bas_share": L["bas"] / (L["bas"] + L["runner"]) if L["bas"] + L["runner"] > 0 else np.nan,
        "anastomoses": float(raw["ana_new"][:, :j + 1].sum()),                       # seit t*
        "front_position": float((r * tips).sum() / tips.sum()) if tips.sum() > 0 else np.nan,   # spitzengewichteter Radius [mm]
    }


def wave_outputs(raw):
    # Fenster bis zum letzten Auswertezeitpunkt
    keep = np.concatenate(([True], raw["times"][1:] * 24 <= max(config.SNAPSHOT_HOURS) + 1e-6))
    sim = {k: (v[:, keep] if isinstance(v, np.ndarray) and v.ndim == 2 else v) for k, v in raw.items()}
    sim["times"] = raw["times"][keep]
    r = ring_centers_mm(raw["location"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fits = amfd.collectArrivalFits([sim])
        C = amfd.encounterC([sim], min_tips=config.ENCOUNTER_MIN_TIPS)["total"]
    rows = []
    for f in fits:
        K1, lam_len, t_n = f["length"]
        K2, lam_tip = (f["tips"][0], f["tips"][1]) if f["tips"] is not None else (np.nan, np.nan)
        rows.append({"r": r[f["ring"]], "t_n": t_n, "K1": K1, "R2_length": f["R2_length"],
                     "saturation": f["rho_data"][-1] / K1, "K2": K2, "lam_tip": lam_tip, "R2_tips": f["R2_tips"]})
    df = pd.DataFrame(rows, columns=["r", "t_n", "K1", "R2_length", "saturation", "K2", "lam_tip", "R2_tips"])
    good = df[df.R2_length >= config.R2_MIN_LENGTH]
    good_K1 = good[good.saturation >= config.SATURATION_MIN]
    good_tips = good[good.R2_tips >= config.R2_MIN_TIPS]
    v = np.nan
    if len(good) >= config.MIN_RINGS and good.r.nunique() >= 2:
        slope = stats.theilslopes(good.t_n, good.r)[0]      # [h/mm]
        v = 24 / slope if slope > 0 else np.nan             # [mm/d]
    med = lambda d, x: float(np.median(x)) if len(d) >= config.MIN_RINGS else np.nan
    return {
        "front_velocity": v,                                            # [mm/d]
        "K1": med(good_K1, good_K1.K1),                                 # [µm mm^-2]
        "K2": med(good_tips, good_tips.K2),                             # [mm^-2]
        "wave_width": med(good_tips, FWHM_BELL / good_tips.lam_tip),    # Dauer der Spitzenwelle im Ring [h]
        "encounter_C": float(C),                                        # f = C rho [mm^2 µm^-1 h^-1]
        "n_wave_rings": len(good), "n_K1_rings": len(good_K1), "n_tip_rings": len(good_tips),
    }


def evaluate_run(run_id, force=False):
    rdir = config.RUNS_DIR / run_id
    meta_path = rdir / "meta.json"
    if not meta_path.exists():
        print(f"{run_id}: noch nicht simuliert"); return
    meta = json.loads(meta_path.read_text())
    out_path = rdir / "outputs.csv"
    if meta.get("sim_status") != "ok":
        print(f"{run_id}: Simulation fehlgeschlagen, keine Auswertung"); return
    if out_path.exists() and not force:
        print(f"{run_id}: bereits ausgewertet"); return
    raw = model_hook.load_raw(rdir / meta["snapshots"][0])
    rows = []
    for h in config.SNAPSHOT_HOURS:
        for k, v in snapshot_outputs(raw, h).items():
            rows.append({"run_id": run_id, "time": h, "region": "all", "output": k, "value": v})
    for k, v in wave_outputs(raw).items():
        rows.append({"run_id": run_id, "time": max(config.SNAPSHOT_HOURS), "region": "all", "output": k, "value": v})
    rows.append({"run_id": run_id, "time": 0, "region": "all", "output": "wait_days", "value": float(raw["wait_days"])})
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f"{run_id}: {len(rows)} Werte")


if __name__ == "__main__":
    design = pd.read_csv(sys.argv[1])
    force = "--force" in sys.argv
    for i in parse_rows(sys.argv[2]):
        if i < len(design):
            evaluate_run(design.iloc[i]["run_id"], force)
