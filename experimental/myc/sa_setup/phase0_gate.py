"""
Automatic check after Phase 0 and the long test, before Morris and the calibration are started (slurm/gate.sbatch).

    python phase0_gate.py            # Phase 0 and long test
    python phase0_gate.py --no-long  # Phase 0 only
      -> results/gate_report.txt     verdict GO / STOP, every check with its numbers (also the body of the e-mail)
      -> results/between_cv.json     robust CV between the seeds of the noise runs for the targets from one colony
                                     (calibrate.py reads it, B: between_cv)
      exit code 0 = GO, 2 = STOP

STOP only for "badly wrong" (the following stages would compute broken or useless runs); everything else is WARN or INFO
and only goes into the report. Thresholds: GATE_* below. Needs the analyses (analyze_phase0.py, analyze_long.py) first.
"""
import argparse
import json
import sys

import numpy as np
import pandas as pd

import config

# ---- thresholds ----
GATE_MIN_OK_NOISE = 0.8         # share of the noise runs (nominal point) that reach t* and are evaluated
GATE_MAX_BAD_RANGE = 0.5        # share of the random points of the factor space without a result (Morris would lose too much)
GATE_WARN_BAD_RANGE = 0.2
GATE_SIM_LIMIT_H = 8.0          # --time of sim.sbatch
GATE_MAX_RUNTIME_FRAC = 0.75    # 90 % quantile of the runtime of the range runs above this share of the limit: timeouts in Morris
GATE_MIN_DEFINED = 0.3          # share of the noise runs in which a wave target is defined (per run; Morris pools 5 replicates)
GATE_WARN_DEFINED = 0.6
GATE_DX_STOP = 0.10             # |relative difference| dx 0.01 vs. 0.005 of a target, CI without 0
GATE_DX_WARN = 0.02             # CI outside ±2 %
GATE_MIN_OK_LONG = 0.7          # share of the long runs (roots off) with a result
GATE_LONG_STOP = 0.25           # |median relative change| 120 h -> LONG_HOURS of a wave target, paired over the seeds
GATE_LONG_WARN = 0.10
WAVE_TARGETS = ["alpha", "beta", "K1", "K2", "front_velocity", "front_tip_speed"]

R = config.RESULTS_DIR


class Report:
    def __init__(self):
        self.lines, self.stops, self.warns = [], 0, 0

    def add(self, level, text):
        self.stops += level == "STOP"
        self.warns += level == "WARN"
        self.lines.append(f"[{level}] {text}")

    def section(self, title):
        self.lines.append(f"\n== {title} ==")


def robust_cv(v):
    v = pd.Series(v, dtype=float).dropna()
    if len(v) < 5 or v.median() == 0:
        return np.nan
    return float(1.4826 * (v - v.median()).abs().median() / abs(v.median()))


def read(name, rep):
    path = R / name
    if not path.exists():
        rep.add("STOP", f"{name} missing (simulation, collection or analysis failed, see slurmOut/sa_collect-*.out)")
        return None
    return pd.read_csv(path)


def targets():
    """(output, time) of the calibration targets; time = end of the window for the wave quantities."""
    import calibrate
    t_end = max(config.SNAPSHOT_HOURS)
    return [(t.output, t.hour if t.hour is not None else t_end, t.replicates) for t in calibrate.TARGETS]


def check_phase0(rep):
    t_end = max(config.SNAPSHOT_HOURS)
    rep.section("Phase 0")
    design, status, out = read("design_phase0.csv", rep), read("status_phase0.csv", rep), read("outputs_phase0.csv", rep)
    if design is None or status is None or out is None:
        return {}
    status = status.merge(design[["run_id", "group"]], on="run_id")
    out = out.merge(design[["run_id", "group", "seed"]], on="run_id")

    # runs with a result
    for g, lim, warn in (("noise", GATE_MIN_OK_NOISE, None), ("range", 1 - GATE_MAX_BAD_RANGE, 1 - GATE_WARN_BAD_RANGE)):
        s = status[status.group == g]
        ok = (s.status == "ok").mean() if len(s) else np.nan
        counts = ", ".join(f"{k} {v}" for k, v in s.status.value_counts().items())
        level = "STOP" if not ok >= lim else ("WARN" if warn is not None and ok < warn else "OK")
        rep.add(level, f"{g}: {ok:.0%} of {len(s)} runs with a result ({counts}); limit {lim:.0%}")
    # runtime of the range runs (edges of the factor space) against the time limit
    t = status[status.group == "range"].sim_seconds.dropna() / 3600
    if len(t):
        q90 = t.quantile(0.9)
        level = "STOP" if q90 > GATE_MAX_RUNTIME_FRAC * GATE_SIM_LIMIT_H else "OK"
        rep.add(level, f"runtime range runs: median {t.median():.2f} h, 90 % {q90:.2f} h, max {t.max():.2f} h "
                       f"(limit sim.sbatch {GATE_SIM_LIMIT_H:g} h; stop above {GATE_MAX_RUNTIME_FRAC:.0%} of it)")
    t = status[status.group == "noise"].sim_seconds.dropna() / 3600
    if len(t):
        rep.add("INFO", f"runtime nominal (noise): median {t.median():.2f} h, max {t.max():.2f} h")

    # wave targets defined per run at the nominal point
    noise = out[(out.group == "noise") & (out.time == t_end) & (out.region == "all")]
    n_noise = int((status[status.group == "noise"].status == "ok").sum())
    for o in WAVE_TARGETS:
        v = noise[noise.output == o].value
        share = v.notna().sum() / n_noise if n_noise else 0.
        level = "STOP" if share < GATE_MIN_DEFINED else ("WARN" if share < GATE_WARN_DEFINED else "OK")
        rep.add(level, f"{o} at {t_end} h defined in {share:.0%} of the noise runs (median {v.median():.4g})")

    # dx convergence of the targets
    dx = read("phase0_dx.csv", rep)
    if dx is not None and len(dx):
        nom = dx[dx.dx == max(config.PHASE0_DX_VALUES)]
        for o, h, _ in targets():
            r = nom[(nom.output == o) & (nom.time == h) & (nom.region == "all")]
            if r.empty or not np.isfinite(r.rel_diff.iloc[0]):
                continue
            r = r.iloc[0]
            sig = r.rel_ci_lo > 0 or r.rel_ci_hi < 0
            if sig and abs(r.rel_diff) > GATE_DX_STOP:
                level = "STOP"
            elif r.rel_ci_lo > GATE_DX_WARN or r.rel_ci_hi < -GATE_DX_WARN:
                level = "WARN"
            else:
                continue
            rep.add(level, f"dx 0.01 vs. {r.ref_dx:g}: {o}@{h:g}h {r.rel_diff:+.3f} [{r.rel_ci_lo:+.3f}, {r.rel_ci_hi:+.3f}]")
        rep.add("INFO", f"dx: checked {len(targets())} targets (only deviations are listed)")
    elif dx is not None:
        rep.add("WARN", "phase0_dx.csv empty: no dx runs evaluated")

    # rings
    w = R / "phase0_warnings.txt"
    if w.exists():
        for line in w.read_text(encoding="utf-8").splitlines():
            if line.startswith("ATTENTION"):
                rep.add("STOP" if "no ring diagnostics" in line else "WARN", line[len("ATTENTION "):])

    # information only: dt, layer, anaL (model definitions / structural assumptions), CRN
    for name, var, ref in (("phase0_dt.csv", "dt", 24 / config.STEPS_PER_DAY), ("phase0_layer.csv", "layer", config.HYPHAE_LAYER)):
        tab = R / name
        if tab.exists():
            d = pd.read_csv(tab)
            d = d[(d.time == t_end) & (d.region == "all") & d.output.isin(WAVE_TARGETS)]
            for _, r in d.iterrows():
                if r[var] == ref and np.isfinite(r.rel_diff):
                    rep.add("INFO", f"{var} {r[var]:g} vs. {r['ref_' + var]:g}: {r.output} {r.rel_diff:+.3f} "
                                    f"[{r.rel_ci_lo:+.3f}, {r.rel_ci_hi:+.3f}]")
    crn = R / "phase0_crn.csv"
    if crn.exists():
        c = pd.read_csv(crn)
        c = c[c.time == t_end]
        if len(c):
            rep.add("INFO", f"CRN: sd ratio same / independent seeds, median over outputs {c.ratio.median():.2f}")

    # between_cv for the targets from one colony (robust CV between the seeds at the nominal point)
    cvs = {}
    for o, h, reps in targets():
        if reps == 1:
            cvs[o] = robust_cv(out[(out.group == "noise") & (out.output == o) & (out.time == h) & (out.region == "all")].value)
    (R / "between_cv.json").write_text(json.dumps({k: (None if not np.isfinite(v) else round(v, 4)) for k, v in cvs.items()},
                                                  indent=1))
    rep.add("INFO", "between_cv (robust CV between the seeds) -> between_cv.json: "
            + ", ".join(f"{k} {v:.3f}" for k, v in cvs.items()))
    return cvs


def check_long(rep):
    rep.section(f"Long test ({config.LONG_HOURS} h)")
    t_end, t_long = max(config.SNAPSHOT_HOURS), config.LONG_HOURS
    design, status, out = read("design_long.csv", rep), read("status_long.csv", rep), read("outputs_long.csv", rep)
    if design is None or status is None or out is None:
        return
    if "group" not in design:
        design["group"] = "long"
    off = design[design.group == "long"]
    s = status[status.run_id.isin(off.run_id)]
    ok = (s.status == "ok").mean() if len(s) else np.nan
    rep.add("STOP" if not ok >= GATE_MIN_OK_LONG else "OK",
            f"long runs (roots off): {ok:.0%} of {len(s)} with a result; limit {GATE_MIN_OK_LONG:.0%}")
    t = s.sim_seconds.dropna() / 3600
    if len(t):
        rep.add("INFO", f"runtime long runs: median {t.median():.2f} h, max {t.max():.2f} h")
    d = out[out.run_id.isin(off.run_id) & (out.region == "all") & out.output.isin(WAVE_TARGETS)]
    for o in WAVE_TARGETS:
        w = d[d.output == o].pivot_table(index="run_id", columns="time", values="value")
        if t_end not in w or t_long not in w:
            rep.add("WARN", f"{o}: no values at {t_end} and {t_long} h")
            continue
        p = w[[t_end, t_long]].dropna()
        p = p[p[t_end] != 0]
        if len(p) < 3:
            rep.add("WARN", f"{o}: only {len(p)} runs with values at {t_end} and {t_long} h, convergence not judged")
            continue
        rel = float(np.median(p[t_long] / p[t_end] - 1))
        level = "STOP" if abs(rel) > GATE_LONG_STOP else ("WARN" if abs(rel) > GATE_LONG_WARN else "OK")
        rep.add(level, f"{o}: median change {t_end} -> {t_long} h {rel:+.3f} (n = {len(p)}); "
                       f"stop above ±{GATE_LONG_STOP:.0%} (the window HOURS_HYPHAE is too short)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-long", action="store_true")
    a = ap.parse_args()
    rep = Report()
    check_phase0(rep)
    if not a.no_long:
        check_long(rep)
    verdict = "STOP" if rep.stops else "GO"
    head = [f"Phase 0 gate: {verdict}  ({rep.stops} STOP, {rep.warns} WARN)",
            "GO: Morris and the calibration are started automatically. STOP: nothing is started; fix and resubmit "
            "(finished runs are reused)."]
    text = "\n".join(head + rep.lines) + "\n"
    (R / "gate_report.txt").write_text(text, encoding="utf-8")
    print(text)
    sys.exit(0 if verdict == "GO" else 2)


if __name__ == "__main__":
    main()
