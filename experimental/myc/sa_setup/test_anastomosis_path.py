"""
Test of the anastomosis exclusion via the path length L (anaMinDist).

Simulates the Petri dish as model_hook.py (nominal XML) for one L and counts the fusions
by relationship of the two hyphae:
  self        hypha with itself
  parent      tip with its parent
  child       tip with one of its own laterals
  sibling     siblings (same parent)
  same_tree   otherwise same hyphalTreeIndex (same infection point)
  other_tree  different trees

Path length diagnostics: for every fusion the path length along the mycelium between the tip (before the fusion)
and the fusion point, as pathDist in MycorrhizalPlant.cpp (positions as arc length at nodes),
plus the length of the fused hypha and its age at fusion. If the path lengths of the
parent/child/sibling fusions are just above L, they are still neighbourhood artefacts (L too tight);
if they are clearly above, a hypha has bent back (fusion legitimate by the rule).

Expectation:
  --L 0     only the own tip segment is excluded, many self/parent fusions
  --L -1    derived L = distTH/sin(theta) + 2*dx (~0.032 cm), self disappears
  --L 100   practically everything related excluded, similar to the old tree exclusion

One run per L, submit in parallel on the cluster (from the home directory), e.g.:
  sbatch --time=6:00:00 --job-name=ana_L0   CPlantBox/experimental/myc/sa_setup/slurm/sa_run.sbatch test_anastomosis_path.py --L 0
  sbatch --time=6:00:00 --job-name=ana_Lm1  CPlantBox/experimental/myc/sa_setup/slurm/sa_run.sbatch test_anastomosis_path.py --L -1
  sbatch --time=6:00:00 --job-name=ana_L100 CPlantBox/experimental/myc/sa_setup/slurm/sa_run.sbatch test_anastomosis_path.py --L 100
Result: output in slurmOut/sa_run-<job>.out, test_ana/L<L>_seed<seed>.json (summary)
and test_ana/L<L>_seed<seed>_fusions.csv (one row per fusion)
"""
import argparse
import csv
import json
import time
from collections import Counter, defaultdict

import numpy as np

import config
import model_hook

MAX_PATH = 5.0  # [cm], path lengths above are not followed (inf)
THRESHOLDS = (0.032, 0.05, 0.1, 0.2)  # [cm], share of fusions with a path length below


def relation(h, m):
    """Relationship between the fused hypha h and partner m."""
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
    """Path length along the mycelium, as chain/pathDist in MycorrhizalPlant.cpp."""

    def __init__(self):
        self.cache = {}  # organ id -> (node coordinates, arc length per node)

    def nodes(self, o):
        i = o.getId()
        if i not in self.cache:
            x = np.array([[v.x, v.y, v.z] for v in o.getNodes()])
            s = np.concatenate([[0.], np.cumsum(np.linalg.norm(np.diff(x, axis=0), axis=1))]) if len(x) > 1 else np.zeros(len(x))
            self.cache[i] = (x, s)
        return self.cache[i]

    def chain(self, o, k, max_len):
        c = {}  # organ id -> (position s, path d up to there)
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
    """Path length, hyphal length and age for the fusion of h into m (last node of h = fusion point)."""
    n = h.getNumberOfNodes()
    k_tip = n - 2                       # tip before the fusion, node n-1 was appended at the fusion
    xh, sh = pl.nodes(h)
    xm, _ = pl.nodes(m)
    dm = np.linalg.norm(xm - xh[n - 1], axis=1)
    if m.getId() == h.getId():
        dm[n - 1] = np.inf                  # do not find the appended node itself
    k_m = int(np.argmin(dm))
    return {"relation": rel, "st_h": int(h.getParameter("subType")), "st_m": int(m.getParameter("subType")),
            "path": pl.dist(h, k_tip, m, k_m), "len_h": float(sh[k_tip]),
            "age_h": 24 * (h.getNodeCT(n - 1) - h.getNodeCT(0)),  # [h]
            "match": float(dm[k_m])}           # Abstand Fusionsknoten - gefundener Partnerknoten, sollte 0 sein


def summarize(records):
    """Quantiles of the path length [cm] and shares below THRESHOLDS per relation and per subtype pair."""
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
    ap.add_argument("--L", type=float, required=True, help="anaMinDist [cm], <0: derived")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--hours", type=int, default=24, help="simulation time after t* [h]")
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
        except Exception as e:  # diagnostics must not cost the result
            failed += 1
            if failed <= 5:
                print(f"fusion_record failed for hypha {h.getId()}: {e!r}", flush=True)

    length = sum(h.getLength() for h in hyphae)
    result = {"L": args.L, "seed": args.seed, "hours": args.hours, "wall_s": round(wall),
              "anaMinDist_used": stats.get("hyphae_anaMinDist"), "n_hyphae": len(hyphae),
              "length_cm": length, "fusions": sum(counts.values()), "relation": dict(counts),
              "relation_by_subtype": dict(sorted(by_subtype.items()))}
    print(result, flush=True)  # print raw first so that nothing is lost

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
        print(f"summarize failed: {e!r}", flush=True)

    def to_json(o):  # convert numpy arrays and numbers (e.g. from stats) to Python types
        return o.tolist() if hasattr(o, "tolist") else str(o)

    text = json.dumps(result, indent=2, default=to_json)
    print(text)
    (out / f"{stem}.json").write_text(text)


if __name__ == "__main__":
    main()
