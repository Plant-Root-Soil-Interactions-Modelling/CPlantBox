"""
Test des Anastomose-Ausschlusses über die Weglänge L (anaMinDist).

Simuliert die Petrischale wie model_hook.py (Nominal-XML) für ein L und zählt die Fusionen
nach Verwandtschaft der beiden Hyphen:
  self        Hyphe mit sich selbst
  parent      Spitze mit ihrem Elter
  child       Spitze mit einem eigenen Lateral
  sibling     Geschwister (gleicher Elter)
  same_tree   sonst gleicher hyphalTreeIndex (gleicher Infektionspunkt)
  other_tree  verschiedene Bäume

Weglängen-Diagnose: für jede Fusion die Weglänge entlang des Myzels zwischen der Spitze (vor der Fusion)
und dem Fusionspunkt, wie pathDist in MycorrhizalPlant.cpp (Positionen als Bogenlänge an Knoten),
dazu die Länge der fusionierten Hyphe und ihr Alter bei der Fusion. Liegen die Weglängen der
parent/child/sibling-Fusionen knapp über L, sind es weiter Nachbarschaftsartefakte (L zu knapp);
liegen sie deutlich darüber, hat sich eine Hyphe zurückgebogen (Fusion nach der Regel legitim).

Erwartung:
  --L 0     nur das eigene Spitzensegment ist ausgeschlossen, viele self/parent-Fusionen
  --L -1    abgeleitetes L = distTH/sin(theta) + 2*dx (~0.032 cm), self verschwindet
  --L 100   praktisch alles Verwandte ausgeschlossen, ähnlich wie der alte Baumausschluss

Ein Lauf pro L, auf dem Cluster parallel einreichen (aus dem Home-Verzeichnis), z. B.:
  sbatch --time=6:00:00 --job-name=ana_L0   CPlantBox/experimental/myc/sa_setup/slurm/sa_run.sbatch test_anastomosis_path.py --L 0
  sbatch --time=6:00:00 --job-name=ana_Lm1  CPlantBox/experimental/myc/sa_setup/slurm/sa_run.sbatch test_anastomosis_path.py --L -1
  sbatch --time=6:00:00 --job-name=ana_L100 CPlantBox/experimental/myc/sa_setup/slurm/sa_run.sbatch test_anastomosis_path.py --L 100
Ergebnis: Ausgabe in slurmOut/sa_run-<job>.out, test_ana/L<L>_seed<seed>.json (Zusammenfassung)
und test_ana/L<L>_seed<seed>_fusions.csv (eine Zeile pro Fusion)
"""
import argparse
import csv
import json
import time
from collections import Counter, defaultdict

import numpy as np

import config
import model_hook

MAX_PATH = 5.0  # [cm], Weglängen darüber werden nicht verfolgt (inf)
THRESHOLDS = (0.032, 0.05, 0.1, 0.2)  # [cm], Anteil der Fusionen mit Weglänge darunter


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


class PathLength:
    """Weglänge entlang des Myzels, wie chain/pathDist in MycorrhizalPlant.cpp."""

    def __init__(self):
        self.cache = {}  # Organ-Id -> (Knotenkoordinaten, Bogenlänge je Knoten)

    def nodes(self, o):
        i = o.getId()
        if i not in self.cache:
            x = np.array([[v.x, v.y, v.z] for v in o.getNodes()])
            s = np.concatenate([[0.], np.cumsum(np.linalg.norm(np.diff(x, axis=0), axis=1))]) if len(x) > 1 else np.zeros(len(x))
            self.cache[i] = (x, s)
        return self.cache[i]

    def chain(self, o, k, max_len):
        c = {}  # Organ-Id -> (Position s, Weg d bis dorthin)
        cur, s, d = o, self.nodes(o)[1][k], 0.
        while cur is not None:
            c[cur.getId()] = (s, d)
            d += s
            if d >= max_len:
                break
            p = cur.getParent()
            if p is None:
                break
            s = self.nodes(p)[1][cur.getParentNI()]
            cur = p
        return c

    def dist(self, a, ka, b, kb, max_len=MAX_PATH):
        ca, cb = self.chain(a, ka, max_len), self.chain(b, kb, max_len)
        best = np.inf
        for i, (sa, da) in ca.items():
            if i in cb:
                sb, db = cb[i]
                best = min(best, da + db + abs(sa - sb))
        return best


def fusion_record(pl, h, m, rel):
    """Weglänge, Hyphenlänge und Alter für die Fusion von h in m (letzter Knoten von h = Fusionspunkt)."""
    n = h.getNumberOfNodes()
    k_tip = n - 2                       # Spitze vor der Fusion, Knoten n-1 wurde bei der Fusion angehängt
    xh, sh = pl.nodes(h)
    xm, _ = pl.nodes(m)
    dm = np.linalg.norm(xm - xh[n - 1], axis=1)
    if m.getId() == h.getId():
        dm[n - 1] = np.inf                  # den angehängten Knoten selbst nicht finden
    k_m = int(np.argmin(dm))
    return {"relation": rel, "st_h": int(h.getParameter("subType")), "st_m": int(m.getParameter("subType")),
            "path": pl.dist(h, k_tip, m, k_m), "len_h": float(sh[k_tip]),
            "age_h": 24 * (h.getNodeCT(n - 1) - h.getNodeCT(0)),  # [h]
            "match": float(dm[k_m])}           # Abstand Fusionsknoten - gefundener Partnerknoten, sollte 0 sein


def summarize(records):
    """Quantile der Weglänge [cm] und Anteile unter THRESHOLDS je Relation und je Subtyp-Paar."""
    groups = defaultdict(list)
    for r in records:
        groups[r["relation"]].append(r)
        groups[f'{r["st_h"]}->{r["st_m"]} {r["relation"]}'].append(r)
    out = {}
    for key, rs in sorted(groups.items()):
        p = np.array([r["path"] for r in rs])
        fin = p[np.isfinite(p)]
        q = np.quantile(fin, [0., 0.1, 0.5, 0.9]) if len(fin) else [np.nan] * 4
        out[key] = {"n": len(rs), "path_min": round(float(q[0]), 4), "path_q10": round(float(q[1]), 4),
                    "path_median": round(float(q[2]), 4), "path_q90": round(float(q[3]), 4),
                    f"path>{MAX_PATH:g}": int(np.sum(~np.isfinite(p))),
                    **{f"frac<{t:g}": round(float(np.mean(p < t)), 3) for t in THRESHOLDS},
                    "len_h_median": round(float(np.median([r["len_h"] for r in rs])), 4),
                    "age_h_median": round(float(np.median([r["age_h"] for r in rs])), 2)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=float, required=True, help="anaMinDist [cm], <0: abgeleitet")
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
    records, failed = [], 0
    pl = PathLength()
    hyphae = mycp.getOrgans(pb.hyphae)
    for h in hyphae:
        m = h.getMergedHyphae()
        if m is None:
            continue
        r = relation(h, m)
        counts[r] += 1
        by_subtype[f"{h.getParameter('subType'):.0f}->{m.getParameter('subType'):.0f} {r}"] += 1
        try:
            records.append(fusion_record(pl, h, m, r))
        except Exception as e:  # Diagnose darf das Ergebnis nicht kosten
            failed += 1
            if failed <= 5:
                print(f"fusion_record fehlgeschlagen für Hyphe {h.getId()}: {e!r}", flush=True)

    length = sum(h.getLength() for h in hyphae)
    result = {"L": args.L, "seed": args.seed, "hours": args.hours, "wall_s": round(wall),
              "anaMinDist_used": stats.get("hyphae_anaMinDist"), "n_hyphae": len(hyphae),
              "length_cm": length, "fusions": sum(counts.values()), "relation": dict(counts),
              "relation_by_subtype": dict(sorted(by_subtype.items()))}
    print(result, flush=True)  # zuerst roh ausgeben, damit nichts verloren geht

    out = config.ROOT / "test_ana"
    out.mkdir(exist_ok=True)
    stem = f"L{args.L:g}_seed{args.seed}"
    if records:
        with open(out / f"{stem}_fusions.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(records[0]))
            w.writeheader()
            w.writerows(records)
    try:
        result["path_diagnostics"] = summarize(records)
        result["path_failed"] = failed
        result["match_max"] = max((r["match"] for r in records), default=None)
    except Exception as e:
        print(f"summarize fehlgeschlagen: {e!r}", flush=True)

    def to_json(o):  # numpy-Arrays und -Zahlen (z. B. aus stats) in Python-Typen umwandeln
        return o.tolist() if hasattr(o, "tolist") else str(o)

    text = json.dumps(result, indent=2, default=to_json)
    print(text)
    (out / f"{stem}.json").write_text(text)


if __name__ == "__main__":
    main()
