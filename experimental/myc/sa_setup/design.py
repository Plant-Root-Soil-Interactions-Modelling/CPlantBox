"""
Erzeugt die Designtabellen (eine Zeile = ein Simulationslauf).

    python design.py phase0     -> results/design_phase0.csv
    python design.py morris     -> results/design_morris.csv (+ morris_X.npy)

Danach: `python run_one.py <design.csv> <zeile>` bzw. per SLURM-Array.
"""
import sys

import numpy as np
import pandas as pd
from SALib.sample import morris as morris_sample

import config
from xml_tools import nominal_values


def morris_problem():
    """SALib-Problemdefinition auf dem normierten Faktorraum [0,1]^k.
    Die Umrechnung in physikalische Werte (inkl. log) macht `to_physical`."""
    k = len(config.FACTORS)
    return {"num_vars": k, "names": [f.name for f in config.FACTORS], "bounds": [[0.0, 1.0]] * k}


def to_physical(U):
    """U: (n, k) in [0,1] -> physikalische Werte gemäß Faktorbereichen."""
    out = np.empty_like(U, dtype=float)
    for j, f in enumerate(config.FACTORS):
        if f.log:
            out[:, j] = np.exp(np.log(f.low) + U[:, j] * (np.log(f.high) - np.log(f.low)))
        else:
            out[:, j] = f.low + U[:, j] * (f.high - f.low)
    return out


def phase0():
    nom = nominal_values()
    rows = []
    for i in range(config.PHASE0_NOISE_SEEDS):
        rows.append({"run_id": f"p0_noise_{i:03d}", "group": "noise",
                     "seed": config.BASE_SEED + i, "dx": np.nan, **nom})
    for dx in config.PHASE0_DX_VALUES:
        for i in range(config.PHASE0_DX_SEEDS):
            # gleiche Seeds über alle dx -> Unterschiede sind nicht nur Seed-Rauschen
            rows.append({"run_id": f"p0_dx{dx:g}_{i:03d}", "group": f"dx={dx:g}",
                         "seed": config.BASE_SEED + i, "dx": dx, **nom})
    df = pd.DataFrame(rows)
    path = config.RESULTS_DIR / "design_phase0.csv"
    df.to_csv(path, index=False)
    print(f"{len(df)} Läufe -> {path}")


def phase1_morris():
    problem = morris_problem()
    k, names = problem["num_vars"], problem["names"]
    # Trajektorienauswahl aus MORRIS_CANDIDATES nach maximaler Spreizung (Ruano et al. 2012,
    # schnelle Variante von Campolongo et al. 2007). Zeilen: Trajektorie für Trajektorie à k+1 Punkte.
    X_flat = morris_sample.sample(problem, config.MORRIS_CANDIDATES, num_levels=config.MORRIS_LEVELS,
                                  optimal_trajectories=config.MORRIS_TRAJECTORIES,
                                  local_optimization=True, seed=config.BASE_SEED)
    X = X_flat.reshape(-1, k + 1, k)
    rows = []
    R = config.MORRIS_REPLICATES
    for t in range(X.shape[0]):
        phys = to_physical(X[t])
        for s in range(k + 1):
            for rep in range(R):
                # Common Random Numbers: gleicher Seed entlang einer Trajektorie
                seed = config.BASE_SEED + 100_000 + t * R + rep
                row = {"run_id": f"m_t{t:02d}_s{s:02d}_r{rep}", "traj": t, "step": s,
                       "rep": rep, "seed": seed, "dx": np.nan}
                row.update(dict(zip(names, phys[s])))
                row.update({f"u_{n}": u for n, u in zip(names, X[t, s])})
                rows.append(row)
    df = pd.DataFrame(rows)
    path = config.RESULTS_DIR / "design_morris.csv"
    df.to_csv(path, index=False)
    np.save(config.RESULTS_DIR / "morris_X.npy", X_flat)
    print(f"{X.shape[0]} Trajektorien × {k + 1} Punkte × {R} Replikate = {len(df)} Läufe -> {path}")


if __name__ == "__main__":
    config.RESULTS_DIR.mkdir(exist_ok=True)
    {"phase0": phase0, "morris": phase1_morris}[sys.argv[1]]()
