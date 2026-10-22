"""
Spielzeugmodell NUR zum Testen der Pipeline. Biologisch bedeutungslos.

Liest dieselbe XML wie die echte Simulation, zieht pro Hyphe Parameter aus
N(value, dev), lässt Spitzen als zufällige Wege wachsen (Puller: laterale
Verzweigung alle ln; Runner: Tip-Splitting mit Rate b) und zählt grob Anastomosen.
`b` der Puller wird bewusst NICHT verwendet (b_prob = 0) -> Kontrollfaktor.
"""
import xml.etree.ElementTree as ET

import numpy as np
from scipy.spatial import cKDTree

import config

DT = 0.5
MAX_TIPS = 300


def _read(xml_path):
    root = ET.parse(xml_path).getroot()
    out = {}
    for hy in root.iter("hyphae"):
        st = int(hy.get("subType"))
        d = {}
        for p in hy.iter("parameter"):
            if p.get("name") == "successor":
                d["succ_p"] = float(p.get("percentage").split(",")[0])
            elif p.get("value") is not None:
                d[p.get("name")] = (float(p.get("value")), float(p.get("dev") or 0))
        out[st] = d
    return out


def _draw(rng, prm, name, lo=1e-9):
    m, s = prm[name]
    return max(lo, rng.normal(m, s)) if s > 0 else m


def _unit(v):
    return v / np.linalg.norm(v)


def _rotate_away(rng, d, angle):
    perp = _unit(np.cross(d, rng.normal(size=3)))
    return _unit(np.cos(angle) * d + np.sin(angle) * perp)


def simulate(xml_path, seed, outdir):
    rng = np.random.default_rng(seed)
    P = _read(xml_path)

    def new_hypha(st, pos, d, t0):
        prm = P[st]
        return dict(st=st, pos=np.array(pos, float), d=_unit(d), t0=t0, L=0.0,
                    v=_draw(rng, prm, "v"), lmax=_draw(rng, prm, "lmax"),
                    hlt=_draw(rng, prm, "hlt"), a=_draw(rng, prm, "a"),
                    ln=_draw(rng, prm, "ln", 1e-3), lb=_draw(rng, prm, "lb"),
                    theta=_draw(rng, prm, "theta"), b=_draw(rng, prm, "b"),
                    next_lat=None, alive=True)

    hy = [new_hypha(1, (0, 0, 0), (0, 0, -1), 0.0)]
    segs, rad, sub, born = [], [], [], []
    paths, t = [], 0.0
    for t_snap in config.SNAPSHOT_TIMES:
        while t < t_snap - 1e-9:
            t += DT
            new = []
            for h in hy:
                if not h["alive"]:
                    continue
                if t - h["t0"] > h["hlt"] or h["L"] >= h["lmax"]:
                    h["alive"] = False
                    continue
                step = min(h["v"] * DT, h["lmax"] - h["L"])
                h["d"] = _unit(h["d"] + 0.3 * rng.normal(size=3))
                p0 = h["pos"].copy()
                h["pos"] = p0 + step * h["d"]
                h["L"] += step
                segs.append((p0, h["pos"].copy())); rad.append(h["a"]); sub.append(h["st"]); born.append(t)
                n_active = sum(x["alive"] for x in hy) + len(new)
                if n_active >= MAX_TIPS:
                    continue
                if h["st"] == 1:          # laterale Verzweigung
                    if h["next_lat"] is None:
                        h["next_lat"] = h["lb"]
                    while h["L"] >= h["next_lat"] and n_active < MAX_TIPS:
                        st = 2 if rng.random() < P[1]["succ_p"] else 1
                        new.append(new_hypha(st, h["pos"], _rotate_away(rng, h["d"], h["theta"]), t))
                        h["next_lat"] += h["ln"]; n_active += 1
                else:                     # Tip-Splitting
                    if rng.random() < 1 - np.exp(-h["b"] / 24 * DT):
                        new.append(new_hypha(2, h["pos"], _rotate_away(rng, h["d"], h["theta"] / 2), t))
                        h["d"] = _rotate_away(rng, h["d"], -h["theta"] / 2)
            hy.extend(new)

        S = np.array(segs)
        tips = np.array([h["pos"] for h in hy if h["alive"]]).reshape(-1, 3)
        # grobe Anastomosen: Spitze trifft Segment-Mitte einer anderen Hyphe
        anas = np.empty((0, 3))
        if len(tips) and len(S):
            mids = S.mean(1)
            dist, _ = cKDTree(mids).query(tips, k=3)
            prm_ana = {1: P[1]["ana"][0], 2: P[2]["ana"][0]}
            dth = np.array([P[h["st"]]["distTH"][0] for h in hy if h["alive"]])
            pana = np.array([prm_ana[h["st"]] for h in hy if h["alive"]])
            hit = (dist[:, -1] < dth * 5) & (rng.random(len(tips)) < pana)
            anas = tips[hit]
        p = outdir / f"snapshot_t{t_snap:g}.npz"
        np.savez_compressed(p, seg_start=S[:, 0], seg_end=S[:, 1], radius=np.array(rad),
                            subtype=np.array(sub), tips=tips, anastomoses=anas)
        paths.append(p)
    return paths
