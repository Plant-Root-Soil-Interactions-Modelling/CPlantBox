"""
Kalibrierung der Hyphenparameter auf Zielwerte als Nullstellensuche, Petrischalen-Setup wie die SA.

    python calibrate.py                    # Kalibrierung (unter SLURM: sbatch calibrate.sbatch)
    python calibrate.py --check            # nur Einstellungen, Bereiche und Startpunkt ausgeben
    python calibrate.py --seeds 8 --max-iter 4

Gesucht ist der Parametervektor x mit F(x) = 0, wobei

    F_i(x) = log(y_i(x) / y*_i) / log(1 + tol_i)

y_i: Output i, gepoolt über N_SEEDS Seeds, y*_i: Zielwert, tol_i: akzeptierte relative Abweichung.
|F_i| <= 1 heißt "innerhalb der Toleranz". Die Parameter werden im normierten Raum u in [0,1]
gesucht (log-uniform wie design.to_physical), damit Schritte zwischen Parametern vergleichbar sind.

Verfahren: Broyden (Quasi-Newton, mehrdimensionales Sekantenverfahren)
  1. F am Startpunkt (Nominalwerte der Template-XML).
  2. Startjacobi-Matrix aus Vorwärtsdifferenzen mit großem Schritt FD_STEP (das Signal muss
     deutlich über dem Seed-Rauschen liegen, kleine Differenzenschritte messen nur Rauschen).
  3. Newton-Schritt J du = -F (kleinste Quadrate, falls #Parameter != #Ziele), begrenzt auf MAX_STEP,
     an den Rand geklemmt. Wird |F| nicht kleiner, wird der Schritt halbiert (max. 3 Versuche).
  4. Jede Auswertung aktualisiert J (Broyden-Update), auch abgelehnte Schritte.
Abbruch, wenn alle |F_i| <= 1, nach MAX_ITER Iterationen oder ohne Fortschritt.

Rauschen: Alle Auswertungen verwenden dieselben Seeds (Common Random Numbers), so sind Differenzen
zwischen Parameterpunkten weniger verrauscht. Der Standardfehler je Output kommt aus einem Bootstrap
über die Seeds und wird mit ausgegeben: Liegt |y - y*| unter ~2 SE, ist die Abweichung im Rauschen.

Ausgabe in calibration/: evals.csv (jede Auswertung), jacobian.csv (zuletzt geschätzte Sensitivitäten
dF_i/du_j), best_params.xml (bester Punkt, sonst wie die Template-XML). Neustartfähig: fertige Läufe
(runs/<punkt>/seed_<s>/meta.json) werden nicht neu simuliert, der Algorithmus läuft deterministisch
bis zur Abbruchstelle nach.
"""
import argparse
import hashlib
import json
import os
import time
import traceback
import warnings
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

import config
import evaluate as ev
import model_hook
from xml_tools import nominal_values, write_xml

amfd = ev.amfd


# --------------------------------------------------------------------------
# Einstellungen
# --------------------------------------------------------------------------
@dataclass
class CalParam:
    factor: str                 # Name aus config.FACTORS (bestimmt subTypes und XML-Parameter)
    low: float | None = None    # Suchbereich, None = Bereich der SA (config.FACTORS)
    high: float | None = None
    log: bool | None = None     # None = wie in config.FACTORS


@dataclass
class Target:
    output: str                 # "beta", "K1", "K2" (gepoolte Wellengrößen) oder ein Name aus config.SNAPSHOT_OUTPUTS
    value: float
    tol_rel: float              # akzeptierte relative Abweichung


# Welche Parameter welches Ziel treiben, lernt das Verfahren (Jacobi-Matrix). Die Reihenfolge ist egal,
# aber jedes Ziel braucht mindestens einen Parameter mit deutlichem Einfluss, sonst ist J singulär.
# Erwartete Hauptwirkungen in den Kommentaren. Suchbereiche teils weiter als in der SA.
PARAMS = [
    CalParam("distTH", 0.005, 0.03, log=True),   # beta (Fangradius der Anastomose, beta ~ distTH^2)
    CalParam("R_ln"),                            # K1 (Lateralabstand der Runner -> Filamentdichte)
    CalParam("B_hlt", 3.0, 12.0),                # bas_share (Lebensdauer der BAS -> ihr Längenanteil)
    CalParam("B_b", 0.3, 2.0, log=True),         # K2 (Tip-Splitting der BAS -> Spitzendichte)
]
TARGETS = [
    Target("beta", 23.0, 0.10),        # [µm/h], a = beta rho n (Oyarte Galvez et al. 2025)
    Target("K1", 1400.0, 0.10),        # [µm mm^-2], Sättigung der Filamentdichte (Fig. 2b)
    Target("bas_share", 0.30, 0.10),   # BAS-Anteil an der Hyphenlänge zum Zeitpunkt HOURS
    Target("K2", 0.15, 0.33),          # [mm^-2], Spitzenpeak 0,1–0,2 (Fig. 2b)
]

HOURS = config.HOURS_HYPHAE   # Simulationsdauer nach t* [h]; Wellen-Fits und Snapshot-Outputs beziehen sich darauf
N_SEEDS = 16                  # Seeds je Auswertung (= Kerne des Jobs, dann eine Runde je Auswertung)
SEED0 = config.BASE_SEED + 300_000
FD_STEP = 0.15                # Schritt der Startdifferenzen im normierten Raum
MAX_STEP = 0.25               # max. Schritt je Iteration (Maximumnorm, normierter Raum)
MAX_ITER = 10
MIN_OK = 0.5                  # Anteil der Seeds, die t* erreichen müssen, sonst gilt der Punkt als ungültig
N_BOOT = 200
CAL_DIR = config.ROOT / "calibration"

FMAP = {f.name: f for f in config.FACTORS}
WAVE = ("beta", "K1", "K2")


# --------------------------------------------------------------------------
# Parameterraum
# --------------------------------------------------------------------------
def spec(p: CalParam):
    f = FMAP[p.factor]
    low = f.low if p.low is None else p.low
    high = f.high if p.high is None else p.high
    log = f.log if p.log is None else p.log
    return low, high, log


def to_phys(u):
    x = []
    for p, uj in zip(PARAMS, u):
        low, high, log = spec(p)
        x.append(float(np.exp(np.log(low) + uj * np.log(high / low)) if log else low + uj * (high - low)))
    return x


def to_unit(x):
    u = []
    for p, xj in zip(PARAMS, x):
        low, high, log = spec(p)
        u.append(np.log(xj / low) / np.log(high / low) if log else (xj - low) / (high - low))
    return np.array(u)


# --------------------------------------------------------------------------
# Simulation und Outputs
# --------------------------------------------------------------------------
def run_seed(args):
    """Ein Lauf (im Worker-Prozess). Gibt den Pfad zu raw.npz zurück oder None."""
    values, seed, rdir, hours = args
    meta_path = rdir / "meta.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
        return str(rdir / meta["snapshots"][0]) if meta["sim_status"] == "ok" else None
    rdir.mkdir(parents=True, exist_ok=True)
    xml_path = write_xml(values, rdir / "params.xml")
    meta = {"seed": seed, "hours": hours, "values": values}
    t0 = time.time()
    try:
        snaps = model_hook.simulate(xml_path, seed, rdir, hours=hours)
        meta.update(sim_status="ok", snapshots=[p.name for p in snaps])
    except Exception:
        meta.update(sim_status="failed", error=traceback.format_exc())
    meta["sim_seconds"] = round(time.time() - t0, 1)
    meta_path.write_text(json.dumps(meta, indent=1))
    return str(rdir / meta["snapshots"][0]) if meta["sim_status"] == "ok" else None


def pooled_outputs(sims, tables):
    """Outputs über alle Seeds gepoolt: beta als Summenverhältnis, K1/K2 als Median über die guten
    Ringe aller Seeds (mind. MIN_RINGS), Snapshot-Outputs als Mittel der Seeds."""
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out["beta"] = amfd.anastomosisBeta(sims, min_tips=config.ENCOUNTER_MIN_TIPS)["total"]
    _, good_K1, good_tips = ev.good_rings(pd.concat(tables, ignore_index=True))
    out["K1"] = float(np.median(good_K1.K1)) if len(good_K1) >= config.MIN_RINGS else np.nan
    out["K2"] = float(np.median(good_tips.K2)) if len(good_tips) >= config.MIN_RINGS else np.nan
    snap = [ev.snapshot_outputs(s, HOURS) for s in sims]
    for t in TARGETS:
        if t.output not in WAVE:
            out[t.output] = float(np.nanmean([s[t.output] for s in snap]))
    return out


class Evaluator:
    def __init__(self, pool, n_seeds, hours):
        self.pool, self.seeds, self.hours = pool, [SEED0 + i for i in range(n_seeds)], hours
        self.records = []

    def __call__(self, u, label):
        x = to_phys(u)
        values = {p.factor: xj for p, xj in zip(PARAMS, x)}
        key = hashlib.md5(json.dumps([round(v, 10) for v in x] + [self.hours]).encode()).hexdigest()[:10]
        pdir = CAL_DIR / "runs" / key
        jobs = [(values, s, pdir / f"seed_{s}", self.hours) for s in self.seeds]
        t0 = time.time()
        paths = [p for p in self.pool.map(run_seed, jobs, chunksize=1) if p]
        sims = [ev.truncate(model_hook.load_raw(Path(p)), self.hours) for p in paths]
        tables = [ev.fit_table(s) for s in sims]

        names = [t.output for t in TARGETS]
        ok = len(sims) >= MIN_OK * len(self.seeds)
        y = pooled_outputs(sims, tables) if ok else dict.fromkeys(names, np.nan)
        # Bootstrap über die Seeds -> Standardfehler
        rng = np.random.default_rng(int(key, 16) % 2 ** 32)
        boot = []
        for _ in range(N_BOOT if ok else 0):
            idx = rng.integers(0, len(sims), len(sims))
            boot.append(pooled_outputs([sims[i] for i in idx], [tables[i] for i in idx]))
        se = {n: float(np.nanstd([b[n] for b in boot])) if boot else np.nan for n in names}

        F = np.array([np.log(y[t.output] / t.value) / np.log(1 + t.tol_rel)
                      if np.isfinite(y[t.output]) and y[t.output] > 0 else np.nan for t in TARGETS])
        rec = {"label": label, "point": key, "n_ok": len(sims), "minutes": round((time.time() - t0) / 60, 1),
               **values, **{f"u_{p.factor}": uj for p, uj in zip(PARAMS, u)},
               **{n: y[n] for n in names}, **{f"se_{n}": se[n] for n in names},
               **{f"F_{n}": Fi for n, Fi in zip(names, F)}, "norm_F": float(np.linalg.norm(F))}
        self.records.append(rec)
        pd.DataFrame(self.records).to_csv(CAL_DIR / "evals.csv", index=False)
        report(rec)
        return F


def report(rec):
    pars = ", ".join(f"{p.factor}={rec[p.factor]:.4g}" for p in PARAMS)
    print(f"\n[{rec['label']}] {pars}  ({rec['n_ok']}/{N_SEEDS} Seeds ok, {rec['minutes']} min)")
    for t in TARGETS:
        n = t.output
        flag = "ok" if abs(rec[f"F_{n}"]) <= 1 else ("im Rauschen" if abs(rec[n] - t.value) <= 2 * rec[f"se_{n}"] else "")
        print(f"   {n:>12}: {rec[n]:10.4g} ± {rec[f'se_{n}']:.2g}  (Ziel {t.value:g} ± {t.tol_rel:.0%})  "
              f"F = {rec[f'F_{n}']:+.2f}  {flag}")
    print(f"   |F| = {rec['norm_F']:.3f}", flush=True)


# --------------------------------------------------------------------------
# Broyden
# --------------------------------------------------------------------------
def broyden_update(J, du, dF):
    if not np.all(np.isfinite(dF)) or du @ du == 0:
        return J
    return J + np.outer(dF - J @ du, du) / (du @ du)


def newton_step(J, F, u, eps=1e-6):
    """J du = -F (kleinste Quadrate). Parameter am Rand, die nach außen wollen, werden festgehalten und der
    Schritt mit den übrigen neu gelöst (sonst würde das Klemmen den Newton-Schritt verfälschen)."""
    free = np.ones(len(u), bool)
    while free.any():
        step = np.zeros(len(u))
        step[free] = -np.linalg.lstsq(J[:, free], F, rcond=None)[0]
        blocked = free & (((u <= eps) & (step < 0)) | ((u >= 1 - eps) & (step > 0)))
        if not blocked.any():
            return step
        free &= ~blocked
    return np.zeros(len(u))


def calibrate(evaluate):
    names = [t.output for t in TARGETS]
    nom = nominal_values()
    u = np.clip(to_unit([nom[p.factor] for p in PARAMS]), 0.02, 0.98)
    F = evaluate(u, "start")
    if not np.all(np.isfinite(F)):
        raise SystemExit("Am Startpunkt sind nicht alle Outputs definiert (Fits gescheitert oder zu viele Läufe ohne t*). "
                         "Mehr Seeds (--seeds) oder längere Simulation (--hours) versuchen.")

    # Startjacobi aus Vorwärtsdifferenzen
    J = np.zeros((len(TARGETS), len(PARAMS)))
    for j, p in enumerate(PARAMS):
        h = FD_STEP if u[j] + FD_STEP <= 1 else -FD_STEP
        uj = u.copy(); uj[j] += h
        Fj = evaluate(uj, f"fd_{p.factor}")
        J[:, j] = np.where(np.isfinite(Fj), (Fj - F) / h, 0.)
    save_jacobian(J)

    best_u, best_F = u, F
    for it in range(1, MAX_ITER + 1):
        if np.all(np.abs(F) <= 1):
            print("\nKonvergiert: alle Ziele innerhalb der Toleranz.")
            break
        step = newton_step(J, F, u)
        scale = min(1., MAX_STEP / max(np.max(np.abs(step)), 1e-12))
        accepted = False
        for attempt in range(3):
            u_new = np.clip(u + scale * step, 0., 1.)
            F_new = evaluate(u_new, f"it{it:02d}" + "abc"[attempt] * (attempt > 0))
            J = broyden_update(J, u_new - u, F_new - F)
            if np.all(np.isfinite(F_new)) and np.linalg.norm(F_new) < np.linalg.norm(F):
                u, F, accepted = u_new, F_new, True
                break
            scale /= 2
        save_jacobian(J)
        if np.linalg.norm(F) < np.linalg.norm(best_F):
            best_u, best_F = u, F
        if not accepted:
            print("\nKein Fortschritt mehr (Schritt dreimal halbiert). Rauschen zu groß, Ziel am Rand des Suchbereichs "
                  "oder mit diesen Parametern nicht erreichbar: jacobian.csv und evals.csv ansehen.")
            break
    else:
        print(f"\nMAX_ITER = {MAX_ITER} erreicht.")

    values = dict(zip([p.factor for p in PARAMS], to_phys(best_u)))
    write_xml(values, CAL_DIR / "best_params.xml")
    print("\nBester Punkt:", ", ".join(f"{k} = {v:.4g}" for k, v in values.items()))
    print("  F =", ", ".join(f"{n} {Fi:+.2f}" for n, Fi in zip(names, best_F)))
    at_edge = [p.factor for p, uj in zip(PARAMS, best_u) if uj <= 1e-6 or uj >= 1 - 1e-6]
    if at_edge:
        print("  Am Rand des Suchbereichs:", ", ".join(at_edge), "-> Bereich erweitern?")
    print(f"  -> {CAL_DIR / 'best_params.xml'}")
    return values


def save_jacobian(J):
    pd.DataFrame(J, index=[t.output for t in TARGETS], columns=[p.factor for p in PARAMS]).to_csv(
        CAL_DIR / "jacobian.csv", float_format="%.4g")


def check():
    nom = nominal_values()
    print(f"{len(PARAMS)} Parameter, {len(TARGETS)} Ziele, {N_SEEDS} Seeds, {HOURS} h nach t*, Template {config.TEMPLATE_XML.name}")
    for p in PARAMS:
        low, high, log = spec(p)
        u0 = to_unit([nom[q.factor] for q in PARAMS])[PARAMS.index(p)]
        print(f"  {p.factor:>8}: [{low:g}, {high:g}]{' log' if log else ''}, nominal {nom[p.factor]:.4g} (u = {u0:.2f})")
    for t in TARGETS:
        print(f"  Ziel {t.output:>10} = {t.value:g} ± {t.tol_rel:.0%}")
    if len(PARAMS) != len(TARGETS):
        print("  Hinweis: #Parameter != #Ziele, Schritte nach kleinsten Quadraten")
    print(f"Auswertungen: 1 + {len(PARAMS)} (Differenzen) + bis {MAX_ITER}×3, je {N_SEEDS} Läufe")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--seeds", type=int, default=N_SEEDS)
    ap.add_argument("--hours", type=float, default=HOURS)
    ap.add_argument("--max-iter", type=int, default=MAX_ITER)
    ap.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count())))
    a = ap.parse_args()
    N_SEEDS, HOURS, MAX_ITER = a.seeds, int(a.hours), a.max_iter
    if a.check:
        check()
    else:
        CAL_DIR.mkdir(exist_ok=True)
        check()
        with Pool(min(a.procs, N_SEEDS)) as pool:
            calibrate(Evaluator(pool, N_SEEDS, HOURS))
