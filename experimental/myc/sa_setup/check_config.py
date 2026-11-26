"""
Checks the configuration before anything goes to the cluster.

    python check_config.py
"""
import xml.etree.ElementTree as ET

import config
from xml_tools import nominal_values


def nominal_values_of(param, subtype):
    root = ET.parse(config.TEMPLATE_XML).getroot()
    hy = next(h for h in root.iter("hyphae") if h.get("subType", "").strip() == str(subtype))
    return float(next(p for p in hy.iter("parameter") if p.get("name") == param).get("value"))

MIN_PER_RUN = 100  # simulation + evaluation [min], 120 h after t* (Phase 0 long test: nominal median 1.6 h), for the estimate


def main():
    nom = nominal_values()
    print(f"Nominal model: {config.TEMPLATE_XML}")
    print(f"{'factor':10s} {'nominal':>10s} {'low':>10s} {'high':>10s}  source")
    open_src, problems = 0, []
    for f in config.FACTORS:
        n = nom[f.name]
        src = f.source or "OPEN (placeholder)"
        open_src += f.source is None
        tags = ("   [control]" if f.role == "control" else "") + ("   [calibrated]" if f.calibrate else "") \
            + (f"   [prior {f.prior[0]:.4g}–{f.prior[1]:.4g}]" if f.prior else "")
        print(f"{f.name:12s} {n:10.4g} {f.low:10.4g} {f.high:10.4g}  {src}{tags}")
        if not (f.low <= n <= f.high):
            problems.append(f"{f.name}: nominal value {n} lies outside [{f.low}, {f.high}]")
        if f.low >= f.high:
            problems.append(f"{f.name}: low >= high")
        if f.log and f.low <= 0:
            problems.append(f"{f.name}: log range needs low > 0")
        if f.param == "ana" and not (0 <= f.low and f.high <= 1):
            problems.append(f"{f.name}: probability must lie in [0, 1]")
    # shared factors: all subTypes should have the same nominal value
    for f in config.FACTORS:
        if len(f.subtypes) > 1 and f.nominal is None and f.param not in ("v_cv", "successor"):
            vals = {st: nominal_values_of(f.param, st) for st in f.subtypes}
            if len(set(vals.values())) > 1:
                problems.append(f"{f.name}: shared factor, but different nominal values {vals}")

    # b_prob must be 0/1 for the factor selection to be right
    root = ET.parse(config.TEMPLATE_XML).getroot()
    for hy in root.iter("hyphae"):
        bp = next((p.get("value") for p in hy.iter("parameter") if p.get("name") == "b_prob"), None)
        print(f"\nsubType {hy.get('subType')} ({hy.get('name')}): b_prob = {bp}", end="")
    print()

    k = len(config.FACTORS)
    n0 = (config.PHASE0_NOISE_SEEDS + len(config.PHASE0_DX_VALUES) * config.PHASE0_DX_SEEDS + config.PHASE0_RANGE_POINTS
          + len(config.PHASE0_DT_STEPS_PER_DAY) * config.PHASE0_DT_SEEDS
          + len(config.PHASE0_LAYER_VALUES) * config.PHASE0_LAYER_SEEDS
          + len(config.PHASE0_ANAL_SCALES) * config.PHASE0_ANAL_SEEDS + 2 * config.PHASE0_CRN_SEEDS)
    dt_cost = sum(spd / config.STEPS_PER_DAY for spd in config.PHASE0_DT_STEPS_PER_DAY) * config.PHASE0_DT_SEEDS
    print(f"dt test: {len(config.PHASE0_DT_STEPS_PER_DAY) * config.PHASE0_DT_SEEDS} runs, runtime ≈ {dt_cost:.0f} nominal runs "
          "(grows with the number of steps -> --time in sim.sbatch)")
    n1 = config.MORRIS_TRAJECTORIES * (k + 1) * config.MORRIS_REPLICATES
    print(f"\nPhase 0: {n0} runs  ≈ {n0 * MIN_PER_RUN / 60:.0f} CPU-h")
    print(f"Morris:  {config.MORRIS_TRAJECTORIES} trajectories × {k + 1} points × "
          f"{config.MORRIS_REPLICATES} replicates = {n1} runs  ≈ {n1 * MIN_PER_RUN / 60:.0f} CPU-h")
    ni = config.INOCULUM_TRAJECTORIES * (len(config.INOCULUM_FACTORS) + 1) * config.INOCULUM_REPLICATES
    print(f"Inoculum Morris: {config.INOCULUM_TRAJECTORIES} trajectories × {len(config.INOCULUM_FACTORS) + 1} points × "
          f"{config.INOCULUM_REPLICATES} replicates = {ni} runs  ≈ {ni * MIN_PER_RUN / 60:.0f} CPU-h")
    for f in config.INOCULUM_FACTORS:
        print(f"  {f.name:18s} [{f.low:g}, {f.high:g}]{' log' if f.log else ''}, nominal {f.nominal:g} ({f.target}): {f.source}")
    print(f"dev mode: {config.DEV_MODE}")
    if open_src:
        print(f"\n{open_src} of {k} factor ranges are still placeholders (±25 %).")
    for p in problems:
        print("ERROR:", p)


if __name__ == "__main__":
    main()
