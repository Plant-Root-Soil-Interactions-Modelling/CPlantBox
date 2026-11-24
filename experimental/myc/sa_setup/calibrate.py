"""
Kalibrierung der Hyphenparameter auf Zielwerte als Nullstellensuche, Petrischalen-Setup wie die SA.

    python calibrate.py --check            # Parameter, Ziele, Startpunkt, Aufwand ausgeben
    python calibrate.py                    # Kalibrierung (unter SLURM: sbatch calibrate.sbatch)
    python calibrate.py --seeds 8 --max-nfev 10

Gesucht ist der Parametervektor x mit F(x) = 0, wobei

    F_i(x) = log(y_i(x) / y*_i) / log(1 + tol_i)

y_i: Output i, gepoolt über N_SEEDS Seeds, y*_i: Zielwert, tol_i: akzeptierte relative Abweichung.
|F_i| <= 1 heißt "innerhalb der Toleranz". Die Parameter werden im normierten Raum u in [0,1]
gesucht (log-uniform wie design.to_physical), damit Schritte zwischen Parametern vergleichbar sind.

Löser: scipy.optimize.least_squares (Trust Region Reflective, Box-Grenzen [0,1]). Bei mindestens so vielen
Zielen wie Parametern ist das Minimum von |F|^2 die Nullstelle; scipy.optimize.root kann keine Grenzen.
Mit weniger Zielen als Parametern ist die Lösung nicht eindeutig: dann bricht das Skript ab
(--allow-underdetermined erzwingt den Lauf). Die Jacobi-Matrix kommt aus Vorwärtsdifferenzen mit großem
Schritt DIFF_STEP (das Signal muss deutlich über dem Seed-Rauschen liegen); alle Differenzenpunkte
werden zusammen auf die Kerne verteilt. Die Suche endet vorzeitig, sobald alle |F_i| <= 1.

Rauschen: Alle Auswertungen verwenden dieselben Seeds (Common Random Numbers), so sind Differenzen
zwischen Parameterpunkten weniger verrauscht. Der Standardfehler je Output kommt aus einem Bootstrap
über die Seeds: Liegt |y - y*| unter ~2 SE, ist die Abweichung im Rauschen.

Ausgabe in calibration/: evals.csv (jede Auswertung), jacobian.csv (Sensitivitäten dF_i/du_j, Spalten
nahe 0 = Parameter von den Zielen nicht bestimmt), best_params.xml (bester Punkt). Neustartfähig: fertige
Läufe (runs/<punkt>/seed_<s>/meta.json) werden nicht neu simuliert, der Löser läuft deterministisch nach.
"""
import argparse
import hashlib
import json
import math
import os
import signal
import time
import traceback
import warnings
from collections import Counter
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

import config
import evaluate as ev
import model_hook
from xml_tools import hyphae_nominal, write_hyphae_xml

amfd = ev.amfd
RUNNER, BAS = config.RUNNER, config.BAS
deg = math.radians


# --------------------------------------------------------------------------
# Parameter
# --------------------------------------------------------------------------
@dataclass
class CalParam:
    name: str
    subtypes: tuple             # subTypes der <hyphae>, auf die der Wert geschrieben wird
    param: str                  # name des <parameter>; "successor" = Anteil des ersten Nachfolgers (Runner: BAS)
    low: float                  # Suchbereich
    high: float
    log: bool = False
    also: tuple = ()            # weitere XML-Parameter derselben subTypes, die denselben Wert bekommen (lb = ln)
    prior: tuple | None = None  # (lo, hi): Literaturbereich als ±2σ-Intervall (log-Parameter: in log); None = schwacher Prior
    source: str = ""            # Herkunft des Priors


# Getrennt nach Runner/BAS wie in der SA (config.FACTORS), dazu getrennt tropismN und a (Radius: sehen im
# Experiment sehr verschieden aus, wirkt auf die Anastomose über den Abstand Spitze–Hyphe); sonst gemeinsam.
# Nicht dabei: b der Runner (wirkungslos bei b_prob = 0), ln und lmax der BAS (keine Laterale, lmax ≈ 3 dx),
# la (fest, sehr klein, Template-XML). lb der Runner ist an ln gekoppelt; lb der BAS bleibt fest.
# Bereiche wie in der SA, wo vorhanden; R_ln, B_b, B_hlt, distTH weiter, tropismN, a, successor geschätzt.
# Priors aus den Quellen in config.FACTORS. B_b ohne Literaturprior: die Messung (Oyarte Galvez 2025) ist schon
# das Ziel branching_rate, sie zählte sonst doppelt.
V_SRC = "180–300 µm/h um 240 µm/h (Oyarte Galvez et al. 2025)"
THETA_SRC = "60° (Schnepf 2016) bis ≈ 70° (Oyarte Galvez 2025), auf 55–75° erweitert (Winkeldefinitionen verschieden)"
PARAMS = [
    # ---- RunnerHyphae (subType 1) ----
    CalParam("R_v", (RUNNER,), "v", 0.432, 0.72, prior=(0.432, 0.72), source=V_SRC),
    CalParam("R_hlt", (RUNNER,), "hlt", 4.0, 10.0, prior=(4.0, 12.5),
             source="SA-Bereich Schnepf 2016 4–10 d, Sterberaten Schnepf 2008 4–12,5 d"),
    CalParam("R_lmax", (RUNNER,), "lmax", 1.5, 5.0),
    CalParam("R_ln", (RUNNER,), "ln", 0.05, 0.6, log=True, also=("lb",), prior=(0.1, 0.6),
             source="~1 mm (Giovannetti, in Schnepf 2008) bis ~6 mm (Oyarte Galvez 2025)"),
    CalParam("R_theta", (RUNNER,), "theta", deg(50), deg(80), prior=(deg(55), deg(75)), source=THETA_SRC),
    CalParam("R_tropismN", (RUNNER,), "tropismN", 1.0, 8.0),
    CalParam("R_succ_BAS", (RUNNER,), "successor", 0.5, 0.95),   # Anteil BAS an den Lateralen der Runner
    # Radius [cm], nominal 0,005; Untergrenze 10 µm (Bago et al. 1998: Durchmesser der Runner > 20 µm), nur Grenze, kein Prior
    CalParam("R_a", (RUNNER,), "a", 0.001, 0.01, log=True),
    # ---- BAS (subType 2) ----
    CalParam("B_v", (BAS,), "v", 0.432, 0.72, prior=(0.432, 0.72), source=V_SRC),
    CalParam("B_b", (BAS,), "b", 0.3, 2.0, log=True),
    CalParam("B_hlt", (BAS,), "hlt", 3.0, 12.0, prior=(5.0, 9.0), source="BAS degenerieren nach ~7 d (Bago 1998; Kokkoris 2026)"),
    CalParam("B_theta", (BAS,), "theta", deg(50), deg(80), prior=(deg(55), deg(75)), source=THETA_SRC),
    CalParam("B_tropismN", (BAS,), "tropismN", 0.5, 4.0),
    CalParam("B_a", (BAS,), "a", 0.00005, 0.005, log=True, prior=(0.0001, 0.0005),   # Radius [cm], nominal 0,0025
             source="Durchmesser 2–10 µm = Radius 1–5 µm (Bago et al. 1998, Mycologia, doi 10.1080/00275514.1998.12026878)"),
    # ---- gemeinsam ----
    CalParam("ana", (RUNNER, BAS), "ana", 0.5, 1.0),
    CalParam("distTH", (RUNNER, BAS), "distTH", 0.005, 0.03, log=True),
]
# Schwacher Prior für Parameter ohne Literatur: Mitte = Nominalwert, σ = WEAK_PRIOR_SIGMA × Suchbereich
# (in log bei log-Parametern). Hält sie beim Nominalwert, solange die Ziele sie nicht brauchen; None = kein Prior.
WEAK_PRIOR_SIGMA = 0.5


# --------------------------------------------------------------------------
# Ziele
# --------------------------------------------------------------------------
@dataclass
class Target:
    output: str                 # siehe OUTPUTS
    value: float
    tol_rel: float              # akzeptierte relative Abweichung
    hour: float | None = None   # nur Snapshot-Outputs: Zeitpunkt [h nach t*]
    log: bool = True            # F = log(y/y*)/log(1+tol); False: F = (y - y*)/(tol y*), für Outputs, die <= 0 werden können

    def residual(self, y):
        if not np.isfinite(y) or (self.log and y <= 0):
            return np.nan
        return np.log(y / self.value) / np.log(1 + self.tol_rel) if self.log else (y - self.value) / (self.tol_rel * self.value)

    @property
    def label(self):
        return self.output if self.hour is None else f"{self.output}@{self.hour:g}h"


# Verfügbare Outputs. Wellengrößen über das Fenster 0..HOURS, über alle Seeds gepoolt:
WAVE_OUTPUTS = {
    "beta": "[µm/h] a = beta rho n mit a, n der Runner-Spitzen (Summenverhältnis)",
    "branching_rate": "[h^-1] neue Runner-Spitzen je wachsender Runner-Spitze und Stunde",
    "encounter_C": "[mm^2 µm^-1 h^-1] f = C rho",
    "K1": "[µm mm^-2] Sättigung der Filamentdichte (Median der gesättigten Ringe)",
    "K2": "[mm^-2] Peak der Runner-Spitzendichte (Median der Ringe mit Peak im Fenster)",
    "wave_width": "[h] Dauer der Runner-Spitzenwelle im Ring (Median)",
    "front_velocity": "[mm/d] Ringradius gegen Ankunft t_n der Filamentdichte (Sigmoid-Fit, Theil-Sen, alle Seeds)",
    "front_speed": "[mm/d] Steigung des spitzengewichteten Radius ab 10 h (Mittel der Seeds)",
    "root_efficiency_rel": "(E_r - E_r,MST) / (E_r,DT - E_r,MST), 2D, am Ende der Simulation (Mittel der Seeds)",
}
# Snapshot-Outputs (Mittel der Seeds) zum Zeitpunkt Target.hour: config.SNAPSHOT_OUTPUTS
# total_length [cm], active_tips, bas_share, anastomoses (seit t*), front_position [mm], spread_r90 [mm],
# max_tip_density [mm^-2]

# Für eine eindeutige Lösung braucht es mindestens so viele Ziele wie Parameter (hier 16), und jeder
# Parameter muss auf mindestens ein Ziel deutlich wirken. Weitere Messwerte als Target(...) ergänzen.
TARGETS = [
    Target("beta", 23.0, 0.10),                    # Oyarte Galvez et al. 2025
    Target("K1", 1400.0, 0.10),                    # Sättigungsdichte, Fig. 2b, gestrichelte Linie
    Target("bas_share", 0.30, 0.10, hour=60),      # BAS-Anteil an der Hyphenlänge
    Target("K2", 0.15, 0.33),                      # Fig. 2b, Spitzenpeak 0,1–0,2
    Target("branching_rate", 0.04, 0.20),          # 4 % pro Stunde
    Target("front_velocity", 280e-3 * 24, 0.10),   # Wellengeschwindigkeit 280 µm/h = 6,72 mm/d
    Target("root_efficiency_rel", 0.5, 0.20, log=False),   # relative Root Efficiency, Ende (HOURS), kann < 0 sein
    # Target("wave_width", ..., 0.20),
    # Target("total_length", ..., 0.15, hour=30),
    # Target("active_tips", ..., 0.15, hour=30),
    # Target("anastomoses", ..., 0.20, hour=60),
]

HOURS = config.HOURS_HYPHAE   # Simulationsdauer nach t* [h], Fenster der Wellengrößen
N_SEEDS = 16                  # Seeds je Parameterpunkt
SEED0 = config.BASE_SEED + 300_000
DIFF_STEP = 0.15              # Differenzenschritt im normierten Raum
MAX_NFEV = 30                 # max. Iterationen des Lösers (least_squares: max_nfev, ohne Differenzen)
MIN_OK = 0.5                  # Anteil der Seeds, die t* erreichen müssen, sonst gilt der Punkt als ungültig
NAN_PENALTY = 10.             # F_i für undefinierte Outputs (Fit gescheitert), damit der Löser weiterläuft
N_BOOT = 200
CAL_DIR = config.ROOT / "calibration"
# Version der Rohdaten (raw.npz) im Schlüssel der Läufe: erhöhen, wenn model_hook neue Größen speichert, die die Ziele
# brauchen; dann werden alte Läufe nicht wiederverwendet. 2: branch_new_runner
RAW_VERSION = 2


# --------------------------------------------------------------------------
# Parameterraum
# --------------------------------------------------------------------------
def spec_nominal(p):
    return hyphae_nominal(p.subtypes[0], p.param)


def prior_spec(p):
    """(Mitte, σ, Art) in der Skala des Parameters (log oder linear), oder None."""
    g = np.log if p.log else (lambda v: v)
    if p.prior is not None:
        lo, hi = p.prior
        return (g(lo) + g(hi)) / 2, (g(hi) - g(lo)) / 4, "Literatur"
    if WEAK_PRIOR_SIGMA:
        return g(spec_nominal(p)), WEAK_PRIOR_SIGMA * (g(p.high) - g(p.low)), "schwach"
    return None


def prior_residuals(u):
    """Zusatzresiduen (g(x) - Mitte) / σ, je Parameter mit Prior: |z| = 1 heißt eine Standardabweichung."""
    z = []
    for p, xj in zip(PARAMS, to_phys(u)):
        s = prior_spec(p)
        if s:
            z.append(((np.log(xj) if p.log else xj) - s[0]) / s[1])
    return np.array(z)


def prior_jacobian():
    """dz/du ist konstant: g(x) ist linear in u."""
    rows = []
    for j, p in enumerate(PARAMS):
        s = prior_spec(p)
        if s:
            row = np.zeros(len(PARAMS))
            row[j] = (np.log(p.high / p.low) if p.log else p.high - p.low) / s[1]
            rows.append(row)
    return np.array(rows).reshape(-1, len(PARAMS))


def to_phys(u):
    return [float(np.exp(np.log(p.low) + uj * np.log(p.high / p.low)) if p.log else p.low + uj * (p.high - p.low))
            for p, uj in zip(PARAMS, u)]


def to_unit(x):
    return np.array([np.log(xj / p.low) / np.log(p.high / p.low) if p.log else (xj - p.low) / (p.high - p.low)
                     for p, xj in zip(PARAMS, x)])


def xml_updates(x):
    """Physikalische Parameterwerte -> {(subType, name): Wert} für write_hyphae_xml."""
    upd = {}
    for p, xj in zip(PARAMS, x):
        for st in p.subtypes:
            for name in (p.param,) + p.also:
                upd[(st, name)] = xj
    return upd


# --------------------------------------------------------------------------
# Simulation und Outputs
# --------------------------------------------------------------------------
def last_line(tb):
    lines = [l for l in str(tb).strip().splitlines() if l.strip()]
    return lines[-1].strip()[:300] if lines else "?"


def run_seed(args):
    """Ein Lauf (im Worker-Prozess). Gibt (Pfad zu raw.npz oder None, Fehlermeldung oder None) zurück."""
    out = lambda meta: ((str(rdir / meta["snapshots"][0]), None) if meta["sim_status"] == "ok"
                        else (None, last_line(meta.get("error", "?"))))
    updates, seed, rdir, hours = args
    meta_path = rdir / "meta.json"
    if meta_path.exists():
        return out(json.loads(meta_path.read_text()))
    xml_path = write_hyphae_xml(updates, rdir / "params.xml")
    meta = {"seed": seed, "hours": hours, "values": {f"{st}:{n}": v for (st, n), v in updates.items()}}
    t0 = time.time()
    try:
        snaps = model_hook.simulate(xml_path, seed, rdir, hours=hours)
        meta.update(sim_status="ok", snapshots=[p.name for p in snaps])
    except Exception:
        meta.update(sim_status="failed", error=traceback.format_exc())
    meta["sim_seconds"] = round(time.time() - t0, 1)
    meta_path.write_text(json.dumps(meta, indent=1))
    return out(meta)


def prepare(sim):
    """Teure Auswertungen je Lauf einmal (Fits, Snapshots), Bootstrap greift darauf zurück.
    Die Wellen-Fits nehmen die Dichte der Runner-Spitzen: im Experiment (Oyarte Galvez et al. 2025, Methods) zählen nur
    wachsende Spitzen (> 40 Pixel bewegt), BAS-Spitzen kaum; Fig. 2/3 zeigen die Dichte der wachsenden RH-Spitzen."""
    hours = sorted({t.hour for t in TARGETS if t.hour is not None})
    sim_rh = dict(sim, tip_dens=sim["active_tips_runner"] / np.asarray(sim["area"]).reshape(-1, 1))
    return {"sim": sim, "table": ev.fit_table(sim_rh), "front_speed": ev.front_speed(sim, HOURS),
            "snap": {h: ev.snapshot_outputs(sim, h) for h in hours}}


def pooled_outputs(preps):
    sims = [p["sim"] for p in preps]
    good, good_K1, good_tips = ev.good_rings(pd.concat([p["table"] for p in preps], ignore_index=True))
    n_min = config.MIN_RINGS
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wave = {
            # a = beta rho n mit a und n der Runner-Spitzen (wachsende Spitzen wie im Experiment), rho aller Hyphen
            "beta": amfd.anastomosisBeta(sims, min_tips=config.ENCOUNTER_MIN_TIPS)["runner"],
            # neue Runner-Spitzen je wachsender Runner-Spitze und Stunde (alpha in b(n) = alpha n)
            "branching_rate": amfd.branchingRate(sims, min_tips=config.ENCOUNTER_MIN_TIPS, suffix="_runner"),
            "encounter_C": amfd.encounterC(sims, min_tips=config.ENCOUNTER_MIN_TIPS)["total"],
            "K1": np.median(good_K1.K1) if len(good_K1) >= n_min else np.nan,
            "K2": np.median(good_tips.K2) if len(good_tips) >= n_min else np.nan,   # Peak der Runner-Spitzendichte
            "wave_width": np.median(ev.FWHM_BELL / good_tips.lam_tip) if len(good_tips) >= n_min else np.nan,
            # Wellenposition aus dem Sigmoid-Fit der Filamentdichte (Fig. 2d inset): Ringradius gegen Ankunftszeit t_n
            "front_velocity": 24 * amfd.theilSlope(good.t_n, good.r)[0] if len(good) >= n_min else np.nan,
            "front_speed": np.nanmean([p["front_speed"] for p in preps]),
            "root_efficiency_rel": np.nanmean([float(p["sim"].get("root_eff_E_rel", np.nan)) for p in preps]),
        }
    out = {}
    for t in TARGETS:
        if t.hour is None:
            out[t.label] = float(wave[t.output])
        else:
            out[t.label] = float(np.nanmean([p["snap"][t.hour][t.output] for p in preps]))
    return out


class Converged(Exception):
    pass


class Residual:
    """F(u) und Jacobi-Matrix für least_squares. Alle Läufe eines Aufrufs (ein Punkt bzw. alle
    Differenzenpunkte) gehen zusammen auf den Pool."""

    def __init__(self, pool, n_seeds, hours):
        self.pool, self.seeds, self.hours = pool, [SEED0 + i for i in range(n_seeds)], hours
        self.records, self.cache = [], {}
        self.best = None    # (|F|, u, F)
        self.errors, self.eff_errors = Counter(), Counter()   # Fehlermeldungen der Läufe bzw. der Root Efficiency
        self.n_runs, self.J, self.status, self.t_start = 0, None, "läuft", time.time()

    def key(self, u):
        return hashlib.md5(json.dumps([round(v, 10) for v in to_phys(u)] + [self.hours, RAW_VERSION]).encode()).hexdigest()[:10]

    def evaluate(self, us, kind):
        """F für mehrere Punkte; bereits berechnete aus dem Cache."""
        todo = [u for u in us if self.key(u) not in self.cache]
        jobs, owner = [], []
        for i, u in enumerate(todo):
            upd = xml_updates(to_phys(u))
            pdir = CAL_DIR / "runs" / self.key(u)
            jobs += [(upd, s, pdir / f"seed_{s}", self.hours) for s in self.seeds]
            owner += [i] * len(self.seeds)
        t0 = time.time()
        results = self.pool.map(run_seed, jobs, chunksize=1)
        minutes = round((time.time() - t0) / 60, 1)
        self.n_runs += len(results)
        self.errors.update(e for _, e in results if e)
        for i, u in enumerate(todo):
            mine = [p for (p, _), o in zip(results, owner) if o == i and p]
            self.cache[self.key(u)] = self.score(u, mine, kind, minutes)
        return [self.cache[self.key(u)] for u in us]

    def score(self, u, paths, kind, minutes):
        labels = [t.label for t in TARGETS]
        preps = [prepare(ev.truncate(model_hook.load_raw(Path(p)), self.hours)) for p in paths]
        self.eff_errors.update(last_line(p["sim"]["root_eff_error"]) for p in preps if "root_eff_error" in p["sim"])
        ok = len(preps) >= MIN_OK * len(self.seeds)
        y = pooled_outputs(preps) if ok else dict.fromkeys(labels, np.nan)
        rng = np.random.default_rng(int(self.key(u), 16) % 2 ** 32)
        boot = [pooled_outputs([preps[i] for i in rng.integers(0, len(preps), len(preps))])
                for _ in range(N_BOOT if ok else 0)]
        se = {n: float(np.nanstd([b[n] for b in boot])) if boot else np.nan for n in labels}
        F = np.array([t.residual(y[t.label]) for t in TARGETS])

        x = to_phys(u)
        z = prior_residuals(u)
        cost = float(np.sum(F ** 2) + np.sum(z ** 2))   # was least_squares minimiert (× 2)
        rec = {"eval": len(self.records), "kind": kind, "point": self.key(u), "n_ok": len(preps), "minutes": minutes,
               **{p.name: xj for p, xj in zip(PARAMS, x)}, **{f"u_{p.name}": uj for p, uj in zip(PARAMS, u)},
               **y, **{f"se_{n}": se[n] for n in labels}, **{f"F_{n}": Fi for n, Fi in zip(labels, F)},
               "norm_F": float(np.linalg.norm(F)), "norm_prior": float(np.linalg.norm(z)), "cost": cost}
        self.records.append(rec)
        pd.DataFrame(self.records).to_csv(CAL_DIR / "evals.csv", index=False)
        report(rec)
        if np.all(np.isfinite(F)) and (self.best is None or cost < self.best[0]):
            self.best = (cost, np.array(u), F)
        write_report(self)
        return F

    def targets(self, u):
        """Residuen der Ziele (undefinierte -> NAN_PENALTY); Konvergenz, wenn alle Ziele in der Toleranz."""
        F = self.evaluate([u], "point")[0]
        if not np.all(np.isfinite(F)):
            if len(self.records) == 1:
                raise SystemExit("Am Startpunkt sind nicht alle Outputs definiert (Fits gescheitert oder zu viele Läufe "
                                 "ohne t*). Mehr Seeds (--seeds) oder längere Simulation (--hours) versuchen.")
            return np.where(np.isfinite(F), F, NAN_PENALTY)
        if np.all(np.abs(F) <= 1):
            self.best = (0., np.array(u), F)   # Ergebnis ist dieser Punkt, auch wenn ein anderer kleinere Kosten hatte
            raise Converged
        return F

    def __call__(self, u):
        """Ziele und Priors zusammen: das ist der Residuenvektor für least_squares."""
        return np.concatenate([self.targets(u), prior_residuals(u)])

    def jac(self, u):
        """Ziele: Vorwärtsdifferenzen, alle Punkte parallel, am oberen Rand rückwärts. Priors: analytisch."""
        F0 = self.targets(u)
        h = np.where(u + DIFF_STEP <= 1, DIFF_STEP, -DIFF_STEP)
        pts = [u + h[j] * np.eye(len(u))[j] for j in range(len(u))]
        Fs = self.evaluate(pts, "diff")
        J = np.column_stack([(np.where(np.isfinite(Fj), Fj, F0) - F0) / h[j] for j, Fj in enumerate(Fs)])
        save_jacobian(J)
        self.J = J
        write_report(self)
        return np.vstack([J, prior_jacobian()])


def report(rec):
    pars = ", ".join(f"{p.name}={rec[p.name]:.3g}" for p in PARAMS)
    print(f"\n[{rec['eval']} {rec['kind']}] {pars}  ({rec['n_ok']}/{N_SEEDS} Seeds ok)")
    for t in TARGETS:
        n = t.label
        flag = "ok" if abs(rec[f"F_{n}"]) <= 1 else ("im Rauschen" if abs(rec[n] - t.value) <= 2 * rec[f"se_{n}"] else "")
        print(f"   {n:>18}: {rec[n]:10.4g} ± {rec[f'se_{n}']:.2g}  (Ziel {t.value:g} ± {t.tol_rel:.0%})  "
              f"F = {rec[f'F_{n}']:+.2f}  {flag}")
    print(f"   |F| = {rec['norm_F']:.3f}, |Prior| = {rec['norm_prior']:.3f}", flush=True)


HINTS = {
    "K1": "zu wenige gesättigte Ringe mit guten Fits -> längere Simulation (--hours 120) oder mehr Seeds",
    "K2": "zu wenige Ringe mit vollständigem Spitzenpeak -> längere Simulation (--hours 120) oder mehr Seeds",
    "wave_width": "wie K2", "front_velocity": "wie K2",
    "root_efficiency_rel": "Root Efficiency nicht berechnet -> Abschnitt Root-Efficiency-Fehler",
    "beta": "zu wenige Ringe mit >= ENCOUNTER_MIN_TIPS Spitzen", "branching_rate": "wie beta",
}


def write_report(r, status=None):
    """calibration/report.txt: Stand, bester Punkt und Fehleranalyse, nach jeder Auswertung neu geschrieben."""
    if status:
        r.status = status
    labels = [t.label for t in TARGETS]
    recs = pd.DataFrame(r.records)
    L = [f"Kalibrierung, Stand {time.strftime('%Y-%m-%d %H:%M')}, Laufzeit {(time.time() - r.t_start) / 3600:.1f} h",
         f"Status: {r.status}"]
    if r.status == "läuft":
        L.append("  (Ist der Job beendet und steht hier noch 'läuft', wurde er hart abgebrochen: Zeitlimit oder Knoten. "
                 "Einfach erneut einreichen, fertige Läufe werden übernommen.)")
    n_pts = int((recs["kind"] == "point").sum()) if len(recs) else 0
    n_prior = sum(1 for p in PARAMS if prior_spec(p))
    n_lit = sum(1 for p in PARAMS if p.prior is not None)
    L += [f"{len(PARAMS)} Parameter, {len(TARGETS)} Ziele, {n_prior} Priors ({n_lit} Literatur)"
          f"{' (UNTERBESTIMMT)' if len(TARGETS) + n_prior < len(PARAMS) else ''}, {len(r.seeds)} Seeds, {r.hours} h nach t*",
          f"Auswertungen {len(recs)} ({n_pts} Iterationspunkte, {len(recs) - n_pts} Differenzenpunkte), Läufe {r.n_runs}", ""]

    if r.best is not None:
        _, u, F = r.best
        b = recs.loc[(recs[[f"u_{p.name}" for p in PARAMS]] - u).abs().max(axis=1).idxmin()]
        L.append(f"Bester Punkt (Auswertung {b['eval']}, |F| Ziele = {b['norm_F']:.2f}, |z| Priors = {b['norm_prior']:.2f}):")
        z = iter(prior_residuals(u))
        for p, uj in zip(PARAMS, u):
            edge = "  <- AM RAND" if uj <= 1e-6 or uj >= 1 - 1e-6 else ""
            s = prior_spec(p)
            pz = f"Prior {s[2]:<9} z = {next(z):+5.2f}" if s else "kein Prior"
            L.append(f"  {p.name:>11} = {b[p.name]:<10.4g} [{p.low:.4g}, {p.high:.4g}], nominal {spec_nominal(p):<8.4g} {pz}{edge}")
        L.append("Ziele am besten Punkt:")
        for t, Fi in zip(TARGETS, F):
            y, se = b[t.label], b[f"se_{t.label}"]
            flag = "ok" if abs(Fi) <= 1 else ("im Rauschen" if abs(y - t.value) <= 2 * se else "daneben")
            L.append(f"  {t.label:>20}: {y:10.4g} ± {se:<8.2g} Ziel {t.value:g} ± {t.tol_rel:.0%}   F = {Fi:+6.2f}  {flag}")
        L.append("Verlauf der Iterationspunkte (Kosten |F|² + |z|²): " + ", ".join(f"{v:.1f}" for v in recs.loc[recs.kind == "point", "cost"]))
        L.append(f"-> best_params.xml" if (CAL_DIR / "best_params.xml").exists() else "")
    else:
        L.append("Noch kein Punkt mit allen Outputs definiert.")

    L += ["", "PRÜFUNGEN"]
    n_fail = sum(r.errors.values())
    L.append(f"Fehlgeschlagene Läufe: {n_fail} von {r.n_runs}" + (" -> häufigste Fehler:" if n_fail else ""))
    L += [f"  {c:4d} x {msg}" for msg, c in r.errors.most_common(5)]
    if len(recs):
        low = recs[recs.n_ok < MIN_OK * len(r.seeds)]
        if len(low):
            L.append(f"  {len(low)} Auswertungen mit weniger als {MIN_OK:.0%} erfolgreichen Seeds (zählen als ungültig)")
        for t in TARGETS:
            nan = int(recs[t.label].isna().sum())
            if nan:
                L.append(f"Undefiniert: {t.label} in {nan} von {len(recs)} Auswertungen ({HINTS.get(t.output, '')})")
        for t in TARGETS:
            se = recs[f"se_{t.label}"].median()
            if np.isfinite(se) and 2 * se > t.tol_rel * abs(t.value):
                L.append(f"Rauschen: {t.label} SE {se:.3g} > halbe Toleranz -> mehr Seeds oder größere Toleranz")
    if r.best is not None:
        zs = iter(prior_residuals(r.best[1]))
        for p in PARAMS:
            s = prior_spec(p)
            zj = next(zs) if s else 0.
            if s and s[2] == "Literatur" and abs(zj) > 2:
                L.append(f"Widerspruch zu Literatur: {p.name} liegt {zj:+.1f} σ vom Prior ({p.source}) -> die Ziele verlangen "
                         f"einen anderen Wert; Ziel, Modellstruktur oder Prior prüfen")
    if r.eff_errors:
        L.append(f"Root-Efficiency-Fehler in {sum(r.eff_errors.values())} Läufen:")
        L += [f"  {c:4d} x {msg}" for msg, c in r.eff_errors.most_common(3)]
    if r.J is not None:
        J = r.J
        col = np.abs(J).max(axis=0)
        dead = [p.name for p, c in zip(PARAMS, col) if c < 0.3]
        if dead:
            L.append("Ohne deutlichen Einfluss auf die Ziele (max |dF/du| < 0,3), nur durch ihren Prior bestimmt: " + ", ".join(dead))
        for t, row in zip(TARGETS, J):
            top = np.argsort(-np.abs(row))[:3]
            weak = "  <- KEIN Parameter wirkt deutlich" if np.abs(row).max() < 0.3 else ""
            L.append(f"  {t.label:>20} wird getrieben von " + ", ".join(f"{PARAMS[k].name} ({row[k]:+.2f})" for k in top) + weak)
        sv = np.linalg.svd(J, compute_uv=False)
        L.append(f"Jacobi der Ziele: Rang {int((sv > 0.05 * sv[0]).sum())} von {min(J.shape)} "
                 f"(Singulärwerte {', '.join(f'{s:.2g}' for s in sv)}); die Priors ergänzen den Rest")
    if not any(l.startswith(("Fehl", "Undef", "Rausch", "Root", "Ohne", "Widerspruch")) and not l.startswith("Fehlgeschlagene Läufe: 0")
               for l in L[L.index("PRÜFUNGEN"):]):
        L.append("Keine Auffälligkeiten.")
    (CAL_DIR / "report.txt").write_text("\n".join(L) + "\n", encoding="utf-8")


def save_jacobian(J):
    pd.DataFrame(J, index=[t.label for t in TARGETS], columns=[p.name for p in PARAMS]).to_csv(
        CAL_DIR / "jacobian.csv", float_format="%.4g")


# --------------------------------------------------------------------------
# Lösen
# --------------------------------------------------------------------------
def calibrate(residual):
    u0 = np.clip(to_unit([spec_nominal(p) for p in PARAMS]), 0.02, 0.98)
    try:
        res = least_squares(residual, u0, jac=residual.jac, bounds=(0., 1.), method="trf",
                            x_scale=1.0, max_nfev=MAX_NFEV, verbose=1)
        F_end = res.fun[:len(TARGETS)]
        status = (f"beendet: Minimum aus Zielen und Priors gefunden, {int(np.sum(np.abs(F_end) <= 1))} von {len(TARGETS)} "
                  f"Zielen in der Toleranz (least_squares: {res.message})")
        print(f"\nleast_squares: {res.message}")
        print("Nicht alle Ziele innerhalb der Toleranz: Rauschen zu groß, Ziel am Rand des Suchbereichs oder mit "
              "diesen Parametern nicht erreichbar -> report.txt, jacobian.csv und evals.csv ansehen.")
    except Converged:
        status = "KONVERGIERT: alle Ziele innerhalb der Toleranz"
        print("\nKonvergiert: alle Ziele innerhalb der Toleranz.")

    if residual.best is None:
        raise SystemExit("Kein Punkt mit allen Outputs definiert, siehe PRÜFUNGEN.")
    _, u, F = residual.best
    x = to_phys(u)
    write_hyphae_xml(xml_updates(x), CAL_DIR / "best_params.xml")
    print("\nBester Punkt:", ", ".join(f"{p.name} = {xj:.4g}" for p, xj in zip(PARAMS, x)))
    print("  F =", ", ".join(f"{t.label} {Fi:+.2f}" for t, Fi in zip(TARGETS, F)))
    at_edge = [p.name for p, uj in zip(PARAMS, u) if uj <= 1e-6 or uj >= 1 - 1e-6]
    if at_edge:
        print("  Am Rand des Suchbereichs:", ", ".join(at_edge), "-> Bereich erweitern?")
    print(f"  -> {CAL_DIR / 'best_params.xml'}")
    write_report(residual, status)
    return x


def check(procs):
    print(f"{len(PARAMS)} Parameter, {len(TARGETS)} Ziele, {N_SEEDS} Seeds, {HOURS} h nach t*, Template {config.TEMPLATE_XML.name}")
    for p in PARAMS:
        nom = spec_nominal(p)
        u0 = to_unit([spec_nominal(q) for q in PARAMS])[PARAMS.index(p)]
        sts = "+".join(config.SUBTYPE_NAMES[s] for s in p.subtypes)
        s = prior_spec(p)
        if s is None:
            pr = "kein Prior"
        else:
            c, w = (np.exp(s[0]), f"×/÷ {np.exp(s[1]):.3g}") if p.log else (s[0], f"± {s[1]:.3g}")
            pr = f"Prior {s[2]}: {c:.4g} {w}" + (f" ({p.source})" if p.source else "")
        print(f"  {p.name:>11} ({sts:>10} {'='.join((p.param,) + p.also)}): [{p.low:.4g}, {p.high:.4g}]"
              f"{' log' if p.log else ''}, nominal {nom:.4g} (u = {u0:.2f}) | {pr}")
    for t in TARGETS:
        print(f"  Ziel {t.label:>18} = {t.value:g} ± {t.tol_rel:.0%}")
    rounds = math.ceil(len(PARAMS) * N_SEEDS / procs)
    print(f"Je Iteration: 1 Punkt + {len(PARAMS)} Differenzenpunkte = {(len(PARAMS) + 1) * N_SEEDS} Läufe, "
          f"mit {procs} Kernen 1 + {rounds} Runden à ~0,75 h")
    n_prior = sum(1 for p in PARAMS if prior_spec(p))
    print(f"Gleichungen: {len(TARGETS)} Ziele + {n_prior} Priors für {len(PARAMS)} Parameter")
    if len(TARGETS) + n_prior < len(PARAMS):
        print(f"\nUNTERBESTIMMT: die Lösung ist nicht eindeutig. Ziele oder Priors ergänzen "
              f"(WEAK_PRIOR_SIGMA) oder --allow-underdetermined.")
        return False
    return True


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--allow-underdetermined", action="store_true")
    ap.add_argument("--seeds", type=int, default=N_SEEDS)
    ap.add_argument("--hours", type=float, default=HOURS)
    ap.add_argument("--max-nfev", type=int, default=MAX_NFEV)
    ap.add_argument("--procs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count())))
    a = ap.parse_args()
    N_SEEDS, HOURS, MAX_NFEV = a.seeds, int(a.hours), a.max_nfev
    determined = check(a.procs)
    if not a.check:
        if not determined and not a.allow_underdetermined:
            raise SystemExit(1)
        CAL_DIR.mkdir(exist_ok=True)
        with Pool(a.procs) as pool:
            r = Residual(pool, N_SEEDS, HOURS)
            def on_term(*_):   # SLURM schickt bei Zeitlimit SIGTERM
                raise SystemExit("SIGTERM erhalten, vermutlich Zeitlimit des Jobs. Erneut einreichen, fertige Läufe werden übernommen.")
            signal.signal(signal.SIGTERM, on_term)
            try:
                calibrate(r)
            except BaseException as e:
                msg = str(e) if isinstance(e, SystemExit) else traceback.format_exc()
                write_report(r, f"ABGEBROCHEN: {msg}")
                raise
