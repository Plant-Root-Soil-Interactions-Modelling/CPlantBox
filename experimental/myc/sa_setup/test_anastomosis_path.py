"""
Test des Anastomose-Ausschlusses über die Weglänge L (anaMinPath).

Simuliert die Petrischale wie model_hook.py (Nominal-XML) für ein L und zählt die Fusionen
nach Verwandtschaft der beiden Hyphen:
  self        Hyphe mit sich selbst
  parent      Spitze mit ihrem Elter
  child       Spitze mit einem eigenen Lateral
  sibling     Geschwister (gleicher Elter)
  same_tree   sonst gleicher hyphalTreeIndex (gleicher Infektionspunkt)
  other_tree  verschiedene Bäume

Erwartung:
  --L 0     nur das eigene Spitzensegment ist ausgeschlossen, viele parent/sibling-Fusionen
  --L -1    abgeleitetes L = distTH/sin(theta) + 2*dx, parent/sibling kaum noch direkt am Verzweigungspunkt
  --L 100   praktisch alles Verwandte ausgeschlossen, ähnlich wie der alte Baumausschluss

Ein Lauf pro L, auf dem Cluster parallel einreichen (aus dem Home-Verzeichnis), z. B.:
  sbatch --time=6:00:00 --job-name=ana_L0   CPlantBox/experimental/myc/sa_setup/slurm/sa_run.sbatch test_anastomosis_path.py --L 0
  sbatch --time=6:00:00 --job-name=ana_Lm1  CPlantBox/experimental/myc/sa_setup/slurm/sa_run.sbatch test_anastomosis_path.py --L -1
  sbatch --time=6:00:00 --job-name=ana_L100 CPlantBox/experimental/myc/sa_setup/slurm/sa_run.sbatch test_anastomosis_path.py --L 100
Ergebnis: Ausgabe in slurmOut/sa_run-<job>.out und test_ana/L<L>_seed<seed>.json
"""
import argparse
import json
import time
from collections import Counter

import config
import model_hook


def relation(h, m):
    """Verwandtschaft zwischen fusionierter Hyphe h und Partner m."""
    if m.getId() == h.getId():
        return "self"
    hp, mp = h.getParent(), m.getParent()
    if hp is not None and hp.getId() == m.getId():
        return "parent"
    if mp is not None and mp.getId() == h.getId():
        return "child"
    if hp is not None and mp is not None and hp.getId() == mp.getId():
        return "sibling"
    if h.getHyphalTreeIndex() == m.getHyphalTreeIndex():
        return "same_tree"
    return "other_tree"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=float, required=True, help="anaMinPath [cm], <0: abgeleitet")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--hours", type=int, default=24, help="Simulationsdauer nach t* [h]")
    args = ap.parse_args()

    model_hook.import_amf()
    import plantbox as pb
    import AMFSimulationSetup as amfs

    xml = config.TEMPLATE_XML
    dish = amfs.defaultPetridish(nRings=config.N_RINGS, hyphae_layer=config.HYPHAE_LAYER, verbose=False)
    dt, nSteps = amfs.setUpSimulationTime(config.ROOT_DAYS, config.STEPS_PER_DAY)
    stats = {}
    t0 = time.perf_counter()
    mycp, time_star, times = amfs.makeSimulation(
        args.seed, str(xml.parent) + "/", xml.stem, dish["height"], dish["petri_dish"], dish["small_hyphae_dish"],
        dish["half_dish"], dt, nSteps, args.hours, "sim", max_wait_days=config.MAX_WAIT_DAYS, stats=stats,
        turn_off_root_growth=True, hyphae_parameters={"anaMinDist": args.L})
    wall = time.perf_counter() - t0
    if time_star is None:
        raise RuntimeError(f"status {stats['status']}")

    counts, by_subtype = Counter(), Counter()
    hyphae = mycp.getOrgans(pb.hyphae)
    for h in hyphae:
        m = h.getMergedHyphae()
        if m is None:
            continue
        r = relation(h, m)
        counts[r] += 1
        by_subtype[f"{h.getParameter('subType'):.0f}->{m.getParameter('subType'):.0f} {r}"] += 1

    length = sum(h.getLength() for h in hyphae)
    result = {"L": args.L, "seed": args.seed, "hours": args.hours, "wall_s": round(wall),
              "anaMinDist_used": stats.get("hyphae_anaMinDist"), "n_hyphae": len(hyphae),
              "length_cm": length, "fusions": sum(counts.values()), "relation": dict(counts),
              "relation_by_subtype": dict(sorted(by_subtype.items()))}
    print(result, flush=True)  # zuerst roh ausgeben, damit nichts verloren geht

    def to_json(o):  # numpy-Arrays und -Zahlen (z. B. aus stats) in Python-Typen umwandeln
        return o.tolist() if hasattr(o, "tolist") else str(o)

    text = json.dumps(result, indent=2, default=to_json)
    print(text)
    out = config.ROOT / "test_ana"
    out.mkdir(exist_ok=True)
    (out / f"L{args.L:g}_seed{args.seed}.json").write_text(text)


if __name__ == "__main__":
    main()
