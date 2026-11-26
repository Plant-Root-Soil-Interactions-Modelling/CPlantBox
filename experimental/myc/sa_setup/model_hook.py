"""
Interface to the simulation: Petri dish as petridishsimauto.py, with the parameter file of the run.

1. simulate(xml_path, seed, outdir) -> list[Path]
   Root phase, waiting for t* and HOURS_HYPHAE hours of hyphae (AMFSimulationSetup.makeSimulation),
   then the quantities per ring and hour as in petridishsimauto.py. Stored in outdir/raw.npz,
   from which evaluate.py computes all outputs (new outputs without re-simulating).
   Without t* after MAX_WAIT_DAYS an error is raised: the run counts as failed.

2. load_raw(path) -> dict
   Reads raw.npz, "times" in days since t* (as AMFIO.load_petridish_results).

While the environment variable SA_DUMMY=1 is set, a toy model is used
(dummy_model.py), so that the pipeline can be tested without the real simulation.
"""
import os
import sys
from pathlib import Path

import numpy as np

import config

USE_DUMMY = os.environ.get("SA_DUMMY") == "1"


class NoTstar(RuntimeError):
    """The run did not reach t* within MAX_WAIT_DAYS; crossed_cm: hyphal length on the fungal side at the stop."""

    def __init__(self, msg, crossed_cm):
        super().__init__(msg)
        self.crossed_cm = crossed_cm


def import_amf():
    """Make plantbox and the AMF* modules from experimental/myc importable."""
    for p in (config.MYC_DIR, config.CPLANTBOX_DIR / "src", config.CPLANTBOX_DIR):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))


def influx_diagnostics(mycp, x0):
    """Supply of hyphae to the fungal side, per hour over the whole simulation (from the final plant, cheap):
    new hyphae emerging from roots, crossings of hyphal segments over x = x0 (boundary of the fungal side, as the rings)
    inwards and outwards (all and per hyphal sub type), new root length. Bins of 1 h from time 0, "diag_edges" in days
    (absolute simulation time, not relative to t*). Tests whether new hyphae stop coming when the roots are turned off."""
    import plantbox as pb
    import AMFDensities as amfd
    T = mycp.getSimTime()
    edges = np.arange(0., T + 1 / 24 + 1e-9, 1 / 24)
    hist = lambda t, w=None: np.histogram(t, bins=edges, weights=w)[0].astype(float)

    # hyphae whose parent is a root (node 0 = branching point, its creation time = emergence)
    ct_emerge = [o.getNodeCT(0) for o in mycp.getOrgans(pb.hyphae)
                 if o.getParent() is not None and o.getParent().organType() == pb.root]
    ana = pb.SegmentAnalyser(mycp)
    nodes = np.array([[n.x, n.y, n.z] for n in ana.nodes]).reshape(-1, 3)
    segs = np.array([[s.x, s.y] for s in ana.segments], dtype=int).reshape(-1, 2)
    ot = np.array(ana.getParameter("organType")).astype(int)
    st = np.array(ana.getParameter("subType")).astype(int)
    ct = np.array(ana.getParameter("creationTime"))
    length = np.array(ana.getParameter("length"))
    out = {"diag_edges": edges, "diag_root_hyphae_new": hist(ct_emerge),
           "diag_root_length_new": hist(ct[ot == pb.root], length[ot == pb.root])}
    xa, xb = nodes[segs[:, 0], 0], nodes[segs[:, 1], 0]     # segment from the older to the newer node
    hy = ot == pb.hyphae
    inward, outward = hy & (xa < x0) & (xb >= x0), hy & (xa >= x0) & (xb < x0)
    out["diag_cross_in"], out["diag_cross_out"] = hist(ct[inward]), hist(ct[outward])
    for s, key in amfd.HYPHAL_SUBTYPES.items():
        out["diag_cross_in_" + key] = hist(ct[inward & (st == s)])
    return out


def frozen_diagnostics(mycp, x0, edges):
    """Hyphae::simulate simulates the children of a hypha only while it is alive: when a hypha reaches its lifetime hlt,
    its whole subtree stops growing (artefact, the children are younger). Per hour (bins `edges` of influx_diagnostics):
    "diag_dead_hyphae" (hyphae that died until then), "diag_frozen_tips" (tips that are still active but frozen, because
    an ancestor died) and "diag_frozen_tips_fungal" (of them with the tip on the fungal side, x >= x0).
    Death time of a hypha = creation time of node 0 + hlt; freeze time of a tip = earliest death time of its ancestors."""
    import plantbox as pb
    hist = lambda t: np.histogram(t, bins=edges)[0].astype(float)
    death, freeze = {}, {}
    dead_t, frozen_t, frozen_fungal_t = [], [], []
    for o in mycp.getOrgans(pb.hyphae):   # parents before children
        oid = o.getId()
        d = o.getNodeCT(0) + o.getParameter("hlt")
        death[oid] = d if not o.isAlive() else np.inf   # a frozen hypha stays alive (its age is not advanced anymore)
        par = o.getParent()
        if par is not None and par.organType() == pb.hyphae and par.getId() in death:
            freeze[oid] = min(freeze.get(par.getId(), np.inf), death[par.getId()])
        else:
            freeze[oid] = np.inf
        if np.isfinite(death[oid]):
            dead_t.append(death[oid])
        if o.isActive() and np.isfinite(freeze[oid]) and freeze[oid] < death[oid]:
            ft = max(freeze[oid], o.getNodeCT(0))   # frozen from its emergence on if the ancestor was already dead
            frozen_t.append(ft)
            if o.getNode(o.getNumberOfNodes() - 1).x >= x0:
                frozen_fungal_t.append(ft)
    cum = lambda t: np.cumsum(hist(t))
    return {"diag_dead_hyphae": cum(dead_t), "diag_frozen_tips": cum(frozen_t), "diag_frozen_tips_fungal": cum(frozen_fungal_t)}


def opening_tip_diagnostics(mycp, barrier_thickness, layer):
    """Active, alive tips at the end inside the barrier opening (-b/2 < x < b/2), in total and below the hyphal layer
    (z < -layer): the opening is 2 mm deep, the layer behind it only `layer`, tips below it run into the step at x = b/2 and
    stay there. They lie outside the rings (not in the ring outputs). Also the active tips on the fungal side (x >= b/2)."""
    import plantbox as pb
    n_open = n_below = n_fungal = 0
    for o in mycp.getOrgans(pb.hyphae):
        if not (o.isActive() and o.isAlive()):
            continue
        p = o.getNode(o.getNumberOfNodes() - 1)
        if abs(p.x) < barrier_thickness / 2:
            n_open += 1
            n_below += p.z < -layer
        elif p.x >= barrier_thickness / 2:
            n_fungal += 1
    return {"diag_opening_tips": float(n_open), "diag_opening_tips_below_layer": float(n_below),
            "diag_fungal_tips_end": float(n_fungal)}


def tip_speed_diagnostics(mycp, x0, time_star, hours, subtype=None, front_quantile=None, center_y=0.):
    """Growth speed of the growing runner tips per hour after t* (Oyarte Galvez et al. 2025, Fig. 2d: "average growth speed
    of hyphae at the front of the colony at each timestep", ≈ 300 µm/h), from the final plant. For every runner and hour
    (t_k, t_k+1]: speed = arc length grown / time between the node before t_k and the last node created until t_k+1 (exact
    node creation times). Growing tips on the fungal side (x >= x0); the front: tips whose distance from the centre of the
    rings is at least the front_quantile quantile of all growing tips at t_k+1. Returns per hour [µm/h] the mean speed
    at the front, of all growing tips, and their numbers; "diag_speed_hours": end of each interval [h after t*]."""
    import plantbox as pb
    subtype = config.RUNNER if subtype is None else subtype
    q = config.FRONT_QUANTILE if front_quantile is None else front_quantile
    n_h = int(hours)
    edges = time_star + np.arange(n_h + 1) / 24.
    speeds = [[] for _ in range(n_h)]        # (speed, radius) per hour
    for o in mycp.getOrgans(pb.hyphae):
        if int(o.getParameter("subType")) != subtype:
            continue
        n = o.getNumberOfNodes()
        if n < 2:
            continue
        cts = np.array([o.getNodeCT(i) for i in range(n)])
        if cts[-1] <= edges[0] or cts[0] > edges[-1]:
            continue
        s = np.array([o.getLength(i) for i in range(n)])          # arc length from the base [cm]
        ks = np.searchsorted(cts, edges, side="right") - 1          # last node created until each edge
        for k in range(n_h):
            a, b = ks[k], ks[k + 1]
            if b <= a or a < 0 or cts[b] <= cts[a]:
                continue                                            # did not grow in this hour
            p = o.getNode(int(b))
            if p.x < x0:
                continue
            v = (s[b] - s[a]) / (cts[b] - cts[a]) * 1.e4 / 24.     # [cm/d] -> [µm/h]
            speeds[k].append((v, np.hypot(p.x, p.y - center_y)))   # distance from the centre of the rings
    front, allv, n_front, n_all, outside = (np.full(n_h, np.nan) for _ in range(5))
    for k, sv in enumerate(speeds):
        if not sv:
            n_front[k] = n_all[k] = 0
            continue
        v, r = np.array(sv).T
        sel = r >= np.quantile(r, q)
        front[k], allv[k], n_front[k], n_all[k] = v[sel].mean(), v.mean(), sel.sum(), len(v)
        outside[k] = np.mean(r > config.RING_RADIUS)        # share of the growing tips beyond the last ring
    return {"diag_speed_hours": np.arange(1, n_h + 1, dtype=float), "diag_front_tip_speed": front,
            "diag_tip_speed_all": allv, "diag_n_front_tips": n_front, "diag_n_growing_tips": n_all,
            "diag_frac_tips_outside_rings": outside}


def _xy_dir(p, q):
    d = np.array([q.x - p.x, q.y - p.y])
    n = np.linalg.norm(d)
    return d / n if n > 1e-12 else None


def branching_angle_diagnostics(mycp, x0, l_parent=None, l_child=None):
    """Branching angles as in Oyarte Galvez et al. 2025, Extended Data Fig. 11b: absolute angle (0–180°) between the initial
    growth direction of a new-born branch and its parent, in the x-y projection (image), on the fungal side (x >= x0).
    Parent direction: tangent over ±l_parent around the branching point (for a tip split: the last l_parent of the parent);
    child direction: from its first node to the node at arc length l_child. By the sub type of the child: RH (runner
    laterals) and BAS (from runners and from BAS tip splits). Experiment: mean 73° (RH), 80° (BAS)."""
    import plantbox as pb
    lp = config.ANGLE_L_PARENT if l_parent is None else l_parent
    lc = config.ANGLE_L_CHILD if l_child is None else l_child
    angles = {config.RUNNER: [], config.BAS: []}
    for o in mycp.getOrgans(pb.hyphae):
        par = o.getParent()
        st = int(o.getParameter("subType"))
        if par is None or par.organType() != pb.hyphae or st not in angles or o.getNumberOfNodes() < 2:
            continue
        p0 = o.getNode(0)
        if p0.x < x0:
            continue
        s_c = np.array([o.getLength(i) for i in range(o.getNumberOfNodes())])
        kc = int(np.searchsorted(s_c, lc - 1e-9))
        if kc >= len(s_c):
            continue                                      # branch shorter than l_child
        dc = _xy_dir(p0, o.getNode(kc))
        npn = par.getNumberOfNodes()
        k = int(o.getParentNI())                          # index of the branching point on the parent
        s_p = np.array([par.getLength(i) for i in range(npn)])
        k0 = int(np.searchsorted(s_p, s_p[k] - lp - 1e-9))
        k1 = min(int(np.searchsorted(s_p, s_p[k] + lp - 1e-9)), npn - 1)
        dp = _xy_dir(par.getNode(k0), par.getNode(k1)) if k1 > k0 else None
        if dc is None or dp is None:
            continue
        angles[st].append(np.degrees(np.arccos(np.clip(dc @ dp, -1., 1.))))
    out = {}
    for st, key in ((config.RUNNER, "rh"), (config.BAS, "bas")):
        a = np.array(angles[st])
        out.update({f"angle_{key}_mean": float(a.mean()) if len(a) else np.nan,
                    f"angle_{key}_median": float(np.median(a)) if len(a) else np.nan,
                    f"angle_{key}_mean_le100": float(a[a <= 100].mean()) if (a <= 100).any() else np.nan,
                    f"angle_{key}_n": float(len(a))})
    return out


def persistence_diagnostics(mycp, x0, max_lag=None, min_length=None):
    """Persistence length of the runners as in Oyarte Galvez et al. 2025, Extended Data Fig. 11a: per hypha the mean cosine
    of the change of the growth direction <cos Δθ(ΔL)> over all starting positions, for lags ΔL up to max_lag (x-y
    projection), fitted with exp(-ΔL / l_p); returns the median l_p over the hyphae [cm] (experiment ≈ 10 cm, a lower
    bound: the decay is barely visible within 1 cm). Runners on the fungal side with a length >= min_length."""
    import plantbox as pb
    L_max = config.PERSISTENCE_MAX_LAG if max_lag is None else max_lag
    L_min = config.PERSISTENCE_MIN_LENGTH if min_length is None else min_length
    lps = []
    for o in mycp.getOrgans(pb.hyphae):
        if int(o.getParameter("subType")) != config.RUNNER or o.getLength(True) < L_min:
            continue
        nodes = np.array([[o.getNode(i).x, o.getNode(i).y] for i in range(o.getNumberOfNodes())])
        if nodes[:, 0].mean() < x0:
            continue
        seg = np.diff(nodes, axis=0)
        ln = np.linalg.norm(seg, axis=1)
        ok = ln > 1e-9
        u, ln = seg[ok] / ln[ok, None], ln[ok]
        if len(u) < 3:
            continue
        dl = np.median(ln)
        K = min(int(L_max / dl), len(u) - 1)
        lags, c = [], []
        for k in range(1, K + 1):
            c.append(np.mean(np.sum(u[:-k] * u[k:], axis=1)))
            lags.append(k * dl)
        lags, c = np.array(lags), np.array(c)
        good = c > 0
        if good.sum() < 2:
            continue
        # least squares of ln<cos> = -ΔL / l_p through the origin
        denom = -np.sum(lags[good] * np.log(c[good]))
        lps.append(np.sum(lags[good] ** 2) / denom if denom > 1e-12 else np.inf)
    lps = np.array(lps)
    return {"persistence_length_median": float(np.median(lps)) if len(lps) else np.nan, "persistence_n": float(len(lps))}


def simulate(xml_path: Path, seed: int, outdir: Path, hours: float | None = None,
             steps_per_day: int | None = None, roots_on: bool = False, sim_args: dict | None = None) -> list[Path]:
    """hours: simulation time after t* [h], default config.HOURS_HYPHAE (for the long test).
    steps_per_day: time steps per day, default config.STEPS_PER_DAY (for the dt convergence test). The evaluation
    stays hourly (times from makeSimulation), only the numerical step changes.
    roots_on: roots keep growing after t* (as in the experiment); default False: root growth is turned off at t*
    (colonization and hyphal emergence stay on in both cases).
    sim_args: factors with target "sim" (config.INOCULUM_FACTORS, SA_GEOMETRY): "root_days" (duration of the root phase [d],
    default config.ROOT_DAYS), "hyphae_layer" (thickness of the hyphal layer [cm], default config.HYPHAE_LAYER) and
    keyword arguments of makeSimulation (seed_pos_x, hyphal_emergence_density)."""
    hours = int(hours) if hours else config.HOURS_HYPHAE
    steps_per_day = int(steps_per_day) if steps_per_day else config.STEPS_PER_DAY
    sim_args = dict(sim_args or {})
    root_days = float(sim_args.pop("root_days", config.ROOT_DAYS))
    layer = float(sim_args.pop("hyphae_layer", config.HYPHAE_LAYER))   # thickness of the hyphal layer [cm] (factor, layer test)
    if USE_DUMMY:
        import dummy_model
        return dummy_model.simulate(xml_path, seed, outdir, hours, roots_on=roots_on, root_days=root_days, sim_args=sim_args,
                                    layer=layer)

    import_amf()
    import AMFAnalysis as amf
    import AMFDensities as amfd
    import AMFSimulationSetup as amfs

    dish = amfs.defaultPetridish(nRings=config.N_RINGS, hyphae_layer=layer, verbose=False)
    dt, nSteps = amfs.setUpSimulationTime(root_days, steps_per_day)
    stats = {}
    mycp, time_star, times = amfs.makeSimulation(
        seed, str(xml_path.parent) + "/", xml_path.stem, dish["height"], dish["petri_dish"], dish["small_hyphae_dish"],
        dish["half_dish"], dt, nSteps, hours, "sim", max_wait_days=config.MAX_WAIT_DAYS, stats=stats,
        min_growing_tips=config.MIN_GROWING_TIPS, stall_fraction=config.STALL_FRACTION, turn_off_root_growth=not roots_on,
        **sim_args)
    if time_star is None:
        # hyphal length (all hyphae) on the fungal side at the stop: the experiment excluded plates with no more than 1 cm
        # (unhealthy AMF); run_one marks such runs "excluded", the others (more, but no t*) "failed"
        ana = amf.getMycSegmentAnalyser(mycp)
        ana.filter("organType", 5)
        crossed = float(ana.getSummed("length", dish["small_hyphae_dish"]))
        raise NoTstar(f"status {stats['status']} (wait_days {stats.get('wait_days', np.nan):.1f}, "
                      f"hyphal length on the fungal side {crossed:.2f} cm)", crossed)
    if config.SAVE_VTP:
        amf.getMycSegmentAnalyser(mycp, stdwrite=True, filename=str(outdir / "final"), step="", timestamp=False)

    # rings as Oyarte Galvez et al. 2025: around the barrier crossing point of THIS run (median of the crossings until t*),
    # N_RINGS of equal area up to RING_RADIUS (Fig. 2c); "dish": the former rings around the centre of the dish
    x0 = -dish["barrier_thickness"] / 2
    if config.RING_CENTER == "crossing":
        center_y, n_cross = amfs.crossingPoint(mycp, x0, t_max=time_star)
    else:
        center_y, n_cross = 0., 0
    rings, location, area = amfs.ringsAround(center_y, config.RING_RADIUS, config.N_RINGS, dish["small_hyphae_dish"],
                                             layer, dish["barrier_thickness"], dish["radius"])
    area = area[:, None]                                                     # [mm^2]
    # as petridishsimauto.py: rings -> rows, time steps -> columns
    mats = amfd.getParasPerRing(["nodeTips", "anastomosis", "length"], times, mycp, rings,
                                subTypes=list(amfd.HYPHAL_SUBTYPES), increments=["anastomosis"])
    active_tips = amfd.getActiveTipsPerRingPerSubtype(times, mycp, rings)
    zeros = np.zeros((len(rings), len(times) - 1))
    raw = {"seed": seed, "time_star": time_star, "location": location, "area": area[:, 0], "times": times[1:],
           "tip_dens": sum(active_tips.values()) / area, "active_tips": sum(active_tips.values()),
           "ana_new": mats["anastomosis_new"], "length_densities": mats["length"] * 10 / area,
           "branch_new": amfd.getBranchingPerRing(times, mycp, rings),
           "branch_new_runner": amfd.getBranchingPerRing(times, mycp, rings, child_subtype=config.RUNNER),
           "wait_days": stats["wait_days"], "steps_per_day": steps_per_day, "roots_on": roots_on,
           "ring_center_y": center_y, "n_crossings_tstar": n_cross, "ring_radius": config.RING_RADIUS}
    try:   # diagnostics only: an error here should not cost the run
        raw.update(influx_diagnostics(mycp, -dish["barrier_thickness"] / 2))
        raw.update(frozen_diagnostics(mycp, -dish["barrier_thickness"] / 2, raw["diag_edges"]))
        raw.update(opening_tip_diagnostics(mycp, dish["barrier_thickness"], layer))
    except Exception:
        import traceback
        raw["diag_error"] = traceback.format_exc()[-2000:]
    try:   # tip speed at the front (Fig. 2d of Oyarte Galvez et al. 2025), calibration target
        raw.update(tip_speed_diagnostics(mycp, -dish["barrier_thickness"] / 2, time_star, hours, center_y=center_y))
    except Exception:
        import traceback
        raw["speed_error"] = traceback.format_exc()[-2000:]
    try:   # branching angles and persistence length (Extended Data Fig. 11), end of the simulation
        raw.update(branching_angle_diagnostics(mycp, -dish["barrier_thickness"] / 2))
        raw.update(persistence_diagnostics(mycp, -dish["barrier_thickness"] / 2))
    except Exception:
        import traceback
        raw["geometry_error"] = traceback.format_exc()[-2000:]
    if config.ROOT_EFFICIENCY:
        # network quantities of the planar image graph (crossings = nodes) as in the experiment, at IMAGE_HOURS after t* and at the
        # end, from the final plant (segments created until then): root efficiency and BAS share / lengths as the experiment
        # classifies the edges (Oyarte Galvez et al. 2025). Keys "<name>_<h>h"; the end also without suffix (as before)
        import time
        import traceback
        hs = sorted({h for h in config.IMAGE_HOURS + config.EXTRA_HOURS if h <= hours} | {hours})
        raw["image_hours"] = np.array(hs, dtype=float)
        for h in hs:
            t_max = time_star + h / 24.
            sfx = f"_{h:g}h"
            t0 = time.time()
            try:   # an error here (including missing networkx) should only cost this output, not the run
                import AMFNetwork as amfn
                # with Ê_g and Ĉ (global_eff): held back from the calibration, for the validation afterwards
                eff = amfn.relativeRootEfficiency(mycp, dish["small_hyphae_dish"], dish["barrier_thickness"], dim=2, t_max=t_max,
                                                  global_eff=True, max_sources=config.GLOBAL_EFF_MAX_SOURCES, seed=seed)
                raw.update({"root_eff_" + k + sfx: v for k, v in eff.items()})
            except Exception:
                raw["root_eff_E_rel" + sfx] = np.nan
                raw["root_eff_error"] = traceback.format_exc()[-2000:]
            raw["root_eff_seconds" + sfx] = time.time() - t0
            t0 = time.time()
            try:
                import AMFNetwork as amfn
                bs = amfn.basShareImage(mycp, -dish["barrier_thickness"] / 2, t_max=t_max)
                raw.update({"bas_image_" + k + sfx: v for k, v in bs.items()})
            except Exception:
                raw["bas_image_share" + sfx] = np.nan
                raw["bas_image_error"] = traceback.format_exc()[-2000:]
            raw["bas_image_seconds" + sfx] = time.time() - t0
        end = f"_{hours:g}h"
        for k in [k for k in raw if k.endswith(end) and k.startswith(("root_eff_", "bas_image_"))]:
            raw[k[:-len(end)]] = raw[k]
    for st, key in amfd.HYPHAL_SUBTYPES.items():
        raw["length_" + key] = mats[f"length_st{st}"] * 10 / area
        raw["active_tips_" + key] = active_tips.get(st, zeros)
        raw["ana_new_" + key] = mats[f"anastomosis_st{st}_new"]
    path = outdir / "raw.npz"
    np.savez_compressed(path, **raw)
    return [path]


def load_raw(path: Path) -> dict:
    with np.load(path) as d:
        raw = {k: d[k] for k in d.files}
    raw["times"] = raw["times"] - raw["time_star"]
    return raw
