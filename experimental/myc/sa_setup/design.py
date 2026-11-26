"""
Creates the design tables (one row = one simulation run).

    python design.py phase0     -> results/design_phase0.csv (noise, dx, sample of the factor space)
    python design.py long       -> results/design_long.csv (nominal point over LONG_HOURS = 180 h, convergence of the wave fits)
    python design.py morris     -> results/design_morris.csv (+ morris_X.npy)
    python design.py inoculum   -> results/design_inoculum.csv (+ inoculum_X.npy), root/inoculum factors

Then: `python run_one.py <design.csv> <row>` or via a SLURM array.
"""
import sys

import numpy as np
import pandas as pd
from SALib.sample import morris as morris_sample

import config
from xml_tools import nominal_values


# Prefix of the run_ids of new designs: the nominal model changed (anaMinDist derived, A11; parameter list v2, B13),
# runs of the first round (p0_*, m_*) must not be reused under the same run_id
P0 = "p0v3"   # v3: HOURS_HYPHAE = 120 (the v2 runs were simulated over 60 h)


def factors_of(tag, design_path=None):
    """Factor list of a Morris design: "morris" (hyphal parameters; v1 for designs of the first round, which have no
    column factor_set) or "inoculum" (root phase, roots, placement)."""
    if tag == "inoculum":
        return config.INOCULUM_FACTORS
    path = design_path or config.RESULTS_DIR / f"design_{tag}.csv"
    try:
        cols = pd.read_csv(path, nrows=1)
        fs = str(cols["factor_set"].iloc[0]) if "factor_set" in cols else "v1"
    except (FileNotFoundError, pd.errors.EmptyDataError):
        fs = config.FACTOR_SET
    return config.factor_set(fs)


def save(df, path):
    df["factor_set"] = config.FACTOR_SET
    df.to_csv(path, index=False)


def morris_problem(factors=None):
    """SALib problem definition on the normalised factor space [0,1]^k.
    The conversion to physical values (incl. log) is done by `to_physical`."""
    factors = factors or config.FACTORS
    k = len(factors)
    return {"num_vars": k, "names": [f.name for f in factors], "bounds": [[0.0, 1.0]] * k}


def to_physical(U, factors=None):
    """U: (n, k) in [0,1] -> physical values according to the factor ranges."""
    factors = factors or config.FACTORS
    out = np.empty_like(U, dtype=float)
    for j, f in enumerate(factors):
        if f.log:
            out[:, j] = np.exp(np.log(f.low) + U[:, j] * (np.log(f.high) - np.log(f.low)))
        else:
            out[:, j] = f.low + U[:, j] * (f.high - f.low)
    return out


def phase0():
    nom = nominal_values()
    rows = []
    for i in range(config.PHASE0_NOISE_SEEDS):
        rows.append({"run_id": f"{P0}_noise_{i:03d}", "group": "noise",
                     "seed": config.BASE_SEED + i, "dx": np.nan, **nom})
    for dx in config.PHASE0_DX_VALUES:
        for i in range(config.PHASE0_DX_SEEDS):
            # same seeds for all dx -> differences are not just seed noise
            rows.append({"run_id": f"{P0}_dx{dx:g}_{i:03d}", "group": f"dx={dx:g}",
                         "seed": config.BASE_SEED + i, "dx": dx, **nom})
    # random points in the whole factor space: how often does a run not reach t*, how long does it take at the edges
    names = [f.name for f in config.FACTORS]
    U = np.random.default_rng(config.BASE_SEED).random((config.PHASE0_RANGE_POINTS, len(names)))
    for i, phys in enumerate(to_physical(U)):
        rows.append({"run_id": f"{P0}_range_{i:03d}", "group": "range",
                     "seed": config.BASE_SEED + 200_000 + i, "dx": np.nan, **dict(zip(names, phys))})
    # dt convergence: same seeds for all dt; the nominal dt (24 steps per day) is simulated again under its own run_id
    for spd in config.PHASE0_DT_STEPS_PER_DAY:
        for i in range(config.PHASE0_DT_SEEDS):
            rows.append({"run_id": f"{P0}_dt{spd}_{i:03d}", "group": f"dt={24 / spd:g}h",
                         "seed": config.BASE_SEED + i, "dx": np.nan, "steps_per_day": spd, **nom})
    # thickness of the hyphal layer (A4)
    for layer in config.PHASE0_LAYER_VALUES:
        for i in range(config.PHASE0_LAYER_SEEDS):
            rows.append({"run_id": f"{P0}_layer{layer:g}_{i:03d}", "group": f"layer={layer:g}",
                         "seed": config.BASE_SEED + i, "dx": np.nan, "hyphae_layer": layer, **nom})
    # minimal path length for anastomosis as a factor on the derived value (A11)
    for s in config.PHASE0_ANAL_SCALES:
        for i in range(config.PHASE0_ANAL_SEEDS):
            # nom also contains the factor ana_L_scale (= 1): it must come first, otherwise it overwrites s
            rows.append({**nom, "run_id": f"{P0}_anaLs{s:g}_{i:03d}", "group": f"anaL={s:g}",
                         "seed": config.BASE_SEED + i, "dx": np.nan, "ana_L_scale": s})
    # common random numbers (C19): the point moved in one factor, with the seeds of the noise runs and with independent seeds
    f = next(f for f in config.FACTORS if f.name == config.PHASE0_CRN_PARAM)
    g = np.log if f.log else (lambda v: v)
    u0 = (g(nom[f.name]) - g(f.low)) / (g(f.high) - g(f.low))
    u1 = u0 + config.PHASE0_CRN_STEP if u0 + config.PHASE0_CRN_STEP <= 1 else u0 - config.PHASE0_CRN_STEP
    moved = {**nom, f.name: float(to_physical(np.array([[u1]]), [f])[0, 0])}
    for i in range(config.PHASE0_CRN_SEEDS):
        rows.append({"run_id": f"{P0}_crnsame_{i:03d}", "group": "crn_same",
                     "seed": config.BASE_SEED + i, "dx": np.nan, **moved})
        rows.append({"run_id": f"{P0}_crnindep_{i:03d}", "group": "crn_indep",
                     "seed": config.BASE_SEED + 600_000 + i, "dx": np.nan, **moved})
    df = pd.DataFrame(rows)
    path = config.RESULTS_DIR / "design_phase0.csv"
    save(df, path)
    print(f"{len(df)} runs -> {path}")


def long_test():
    """Nominal point over LONG_HOURS h: does the wave develop with more time (wave fits defined often enough)?
    Plus the same seeds with roots growing on after t* (group "long_roots_on"): when do new hyphae stop coming with
    the roots turned off, and from when on do the outputs differ?"""
    nom = nominal_values()
    rows = [{"run_id": f"{P0}_long_{i:03d}", "group": "long", "seed": config.BASE_SEED + i, "dx": np.nan,
             "hours": config.LONG_HOURS, "roots_on": False, **nom} for i in range(config.LONG_SEEDS)]
    rows += [{"run_id": f"{P0}_longroots_{i:03d}", "group": "long_roots_on", "seed": config.BASE_SEED + i, "dx": np.nan,
              "hours": config.LONG_HOURS, "roots_on": True, **nom} for i in range(config.LONG_SEEDS)]
    df = pd.DataFrame(rows)
    path = config.RESULTS_DIR / "design_long.csv"
    save(df, path)
    print(f"{len(df)} runs over {config.LONG_HOURS} h -> {path}")


def morris_design(tag, trajectories, R, seed_offset, prefix):
    factors = config.INOCULUM_FACTORS if tag == "inoculum" else config.FACTORS   # new designs: current list
    problem = morris_problem(factors)
    k, names = problem["num_vars"], problem["names"]
    # trajectory selection from MORRIS_CANDIDATES by maximum spread (Ruano et al. 2012,
    # fast variant of Campolongo et al. 2007). Rows: trajectory by trajectory, k+1 points each.
    X_flat = morris_sample.sample(problem, config.MORRIS_CANDIDATES, num_levels=config.MORRIS_LEVELS,
                                  optimal_trajectories=trajectories,
                                  local_optimization=True, seed=config.BASE_SEED)
    X = X_flat.reshape(-1, k + 1, k)
    # inoculum design: hyphal parameters nominal (run_one needs all of them)
    nom = nominal_values() if tag != "morris" else {}
    rows = []
    for t in range(X.shape[0]):
        phys = to_physical(X[t], factors)
        for s in range(k + 1):
            for rep in range(R):
                # common random numbers: same seed along a trajectory
                seed = config.BASE_SEED + seed_offset + t * R + rep
                row = {"run_id": f"{prefix}_t{t:02d}_s{s:02d}_r{rep}", "traj": t, "step": s,
                       "rep": rep, "seed": seed, "dx": np.nan, **nom}
                row.update(dict(zip(names, phys[s])))
                row.update({f"u_{n}": u for n, u in zip(names, X[t, s])})
                rows.append(row)
    df = pd.DataFrame(rows)
    path = config.RESULTS_DIR / f"design_{tag}.csv"
    save(df, path)
    np.save(config.RESULTS_DIR / f"{tag}_X.npy", X_flat)
    print(f"{X.shape[0]} trajectories × {k + 1} points × {R} replicates = {len(df)} runs -> {path}")


def phase1_morris():
    # second round (parameter list v2): own run_ids and seeds, the first round used "m" and 100_000
    morris_design("morris", config.MORRIS_TRAJECTORIES, config.MORRIS_REPLICATES, 500_000, "m2")


def inoculum_morris():
    """Root phase, root architecture, placement and emergence density, hyphal parameters nominal (config.INOCULUM_FACTORS)."""
    # "inoc3": HOURS_HYPHAE = 120 (the "inoc" runs were simulated over 60 h)
    morris_design("inoculum", config.INOCULUM_TRAJECTORIES, config.INOCULUM_REPLICATES, 400_000, "inoc3")


if __name__ == "__main__":
    config.RESULTS_DIR.mkdir(exist_ok=True)
    {"phase0": phase0, "long": long_test, "morris": phase1_morris, "inoculum": inoculum_morris}[sys.argv[1]]()
