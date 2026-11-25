import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import plantbox as pb
from scipy.optimize import curve_fit
from scipy.spatial import cKDTree

import AMFAnalysis as amf

# Plant -> amounts per ring and time interval (length, tips, anastomoses) -> arrival time fits of the densities

HYPHAL_SUBTYPES = {1: "runner", 2: "bas"} # hyphal sub types of TwoHyphaePlusBAS.xml (RunnerHyphae, BAS)

def getRingIndex(ana, rings):
    # ring index of each segment (-1 = in no ring), uses the segment midpoint like ana.getSummed(name, ring)
    nodes = np.array([[n.x, n.y, n.z] for n in ana.nodes])
    segs = np.array([[s.x, s.y] for s in ana.segments], dtype=int)
    mids = 0.5 * (nodes[segs[:, 0]] + nodes[segs[:, 1]])
    ringIdx = np.full(len(mids), -1)
    for i, m in enumerate(mids):
        p = pb.Vector3d(m[0], m[1], m[2])
        for k, ring in enumerate(rings):
            if ring.getDist(p) < 0:
                ringIdx[i] = k
                break
    return ringIdx

def getParasPerRing(parameters, times, plant, rings, subTypes = None, increments = ()):
    # matrices (rings x len(times)-1) for all parameters at once (same layout as AMFLegacy.getParaSumperRing_old):
    # one SegmentAnalyser, ring of each segment computed once, time intervals binned with numpy
    # subTypes: hyphal sub types (e.g. [1, 2]), additionally returns parameter + "_st" + subType with only the hyphae of that sub type
    # increments: parameters for which additionally parameter + "_new" (and parameter + "_st" + subType + "_new") is returned:
    #             the amount created in each interval (times[i], times[i+1]] instead of the cumulative amount
    ana = amf.getMycSegmentAnalyser(plant)
    ringIdx = getRingIndex(ana, rings)
    ct = np.array(ana.getParameter("creationTime"))
    n = len(times) - 1
    timeIdx = np.searchsorted(times, ct, side="left") # 0: ct <= times[0], i: times[i-1] < ct <= times[i]
    valid = (ringIdx >= 0) & (ct >= 0) & (timeIdx <= n)
    groups = [("", valid)]
    if subTypes is not None:
        ot = np.array(ana.getParameter("organType")).astype(int)
        st = np.array(ana.getParameter("subType")).astype(int)
        groups += [("_st" + str(s), valid & (ot == pb.hyphae) & (st == s)) for s in subTypes]
    mats = {}
    for parameter in parameters:
        values = np.array(ana.getParameter(parameter))
        for suffix, mask in groups:
            summed = np.zeros((len(rings), n + 1))
            np.add.at(summed, (ringIdx[mask], timeIdx[mask]), values[mask])
            if parameter == "nodeTips": # we want to look only at actively growing tips. (Oyarte Galvez 2025 page 3 when discussing Fig 2b)
                mat = summed[:, 1:] # created in (times[i], times[i+1]]
            else:
                mat = np.cumsum(summed, axis=1)[:, 1:] # created until times[i+1]
                mat[:, 0] = summed[:, 1] # first column only counts (times[0], times[1]], as in AMFLegacy.getParaSumperRing_old
            mats[parameter + suffix] = mat
            if parameter in increments:
                mats[parameter + suffix + "_new"] = summed[:, 1:]
    return mats

def getActiveTipsPerRing(times, plant, rings, ot = 5):
    # all sub types together, see getActiveTipsPerRingPerSubtype
    return sum(getActiveTipsPerRingPerSubtype(times, plant, rings, ot).values())

def getActiveTipsPerRingPerSubtype(times, plant, rings, ot = 5):
    # {sub type: matrix} with the matrices of getActiveTipsPerRing for each sub type
    # number of actively growing tips per ring at times[1:] (rings x len(times)-1, same layout as getParasPerRing),
    # reconstructed from the final plant: an organ is growing at time t if it was created (node 0) before t and
    # its last node was created after the start of the interval (times[i], times[i+1] = t]. Its tip at time t is its last node created until t.
    # o.isActive() is not used: a hypha whose ancestor died is not simulated anymore (Hyphae::simulate), it stays active but does not grow
    # (node 1 is not used as start, it gets a new creation time as long as the first segment is shorter than dx and the node is moved)
    # ("nodeTips" in getParasPerRing instead counts organs whose last node was created in the interval, i.e. that stopped growing then)
    t_eval = np.asarray(times)[1:]
    t_prev = np.asarray(times)[:-1]
    mats = {}
    ringOfNode = {} # global node id -> ring index (-1 = in no ring)
    for o in plant.getOrgans(ot):
        st = int(o.getParameter("subType"))
        if st not in mats:
            mats[st] = np.zeros((len(rings), len(t_eval)))
        mat = mats[st]
        n = o.getNumberOfNodes()
        cts = np.array([o.getNodeCT(i) for i in range(n)])
        t_start, t_stop = cts[0], cts[-1] # node 0 is the branching point on the parent, its creation time is the emergence of the organ
        active = (t_eval >= t_start) & (t_prev < t_stop)
        if not active.any():
            continue
        tipIdx = np.searchsorted(cts, t_eval[active], side="right") - 1 # last node created until t
        for j, k in zip(np.nonzero(active)[0], tipIdx):
            nodeId = o.getNodeId(int(k))
            if nodeId not in ringOfNode:
                ringOfNode[nodeId] = -1
                p = o.getNode(int(k))
                for r, ring in enumerate(rings):
                    if ring.getDist(p) < 0:
                        ringOfNode[nodeId] = r
                        break
            r = ringOfNode[nodeId]
            if r >= 0:
                mat[r, j] += 1
    if not mats:
        mats[0] = np.zeros((len(rings), len(t_eval)))
    return mats

def getBranchingPerRing(times, plant, rings, ot = 5, child_subtype = None):
    # branching events per ring in each interval (times[i], times[i+1]] (rings x len(times)-1, same layout as getParasPerRing), at the
    # branching point (node 0 of the child). A tip split creates two children at the same time and place (Hyphae::simulate, b_prob > 0.5),
    # a lateral branch one: children are grouped by (parent, creation time, position), each group is one event.
    # Hyphae emerging from roots are no branching events. child_subtype: only events whose children have this sub type
    # (e.g. 1: new runner tips, as the newborn growing tips of Oyarte Galvez et al. 2025, where BAS tips hardly count as growing)
    times = np.asarray(times)
    mat = np.zeros((len(rings), len(times) - 1))
    events = {}
    for o in plant.getOrgans(ot):
        parent = o.getParent()
        if parent is None or parent.organType() != ot:
            continue
        if child_subtype is not None and int(o.getParameter("subType")) != child_subtype:
            continue
        p = o.getNode(0)
        events[(parent.getId(), round(o.getNodeCT(0), 9), round(p.x, 6), round(p.y, 6), round(p.z, 6))] = p
    for (_, ct, *_), p in events.items():
        j = np.searchsorted(times, ct, side="left") # as getParasPerRing: i means times[i-1] < ct <= times[i]
        if j < 1 or j >= len(times):
            continue
        for r, ring in enumerate(rings):
            if ring.getDist(p) < 0:
                mat[r, j - 1] += 1
                break
    return mat

def branchingRate(simulations, min_tips = 5, suffix = ""):
    # beta: branching events per growing tip and hour [h^-1] of b = beta n (Oyarte Galvez et al. 2025: about 4 %/h; they call it alpha,
    # here beta = branching and alpha = anastomosis, see anastomosisAlpha), as anastomosisAlpha: sum of the events
    # over sum of tips x interval length, over all rings and intervals (first column left out) with at least min_tips growing tips.
    # suffix "_runner": new runner tips per growing runner tip ("branch_new_runner", "active_tips_runner"); nan if the data are missing
    ev_sum, tip_h = 0., 0.
    for sim in simulations:
        if "branch_new" + suffix not in sim:
            return np.nan
        dt_h = np.diff(np.concatenate([[0.], sim["times"]]))[1:] * 24
        tips, br = sim["active_tips" + suffix][:, 1:], sim["branch_new" + suffix][:, 1:]
        ok = tips >= min_tips
        ev_sum += br[ok].sum()
        tip_h += (tips * dt_h)[ok].sum()
    return ev_sum / tip_h if tip_h > 0 else np.nan

def _pointSegmentDistance(p, a, b):
    # distance of the points p to the segments [a, b] (all (n, 3), row by row)
    v = b - a
    vv = np.einsum("ij,ij->i", v, v)
    s = np.clip(np.einsum("ij,ij->i", p - a, v) / np.where(vv > 0, vv, 1.), 0., 1.)
    return np.linalg.norm(p - (a + s[:, None] * v), axis=1)

def _chain(o, k, maxLen):
    # path from node k of organ o down to the base and on through the parents, as chain() in MycorrhizalPlant.cpp:
    # (organ id, position s along that organ [cm], path length d before reaching that organ [cm]), until d > maxLen
    c, cur, s, d = [], o, o.getLength(int(k)), 0.
    while cur is not None:
        c.append((cur.getId(), s, d))
        d += s
        if d > maxLen:
            break
        par = cur.getParent()
        if par is None:
            break
        s, cur = par.getLength(int(cur.getParentNI())), par
    return c

def _pathDist(chainA, chainB):
    # path length along the mycelium between the two points of the chains, as pathDist() in MycorrhizalPlant.cpp (1e100: no common organ)
    best = 1e100
    posB = {}
    for oid, s, d in chainB:
        posB.setdefault(oid, []).append((s, d))
    for oid, sA, dA in chainA:
        for sB, dB in posB.get(oid, ()):
            best = min(best, dA + dB + abs(sA - sB))
    return best

def anastomosisMinPath(plant, subtype, ot = 5):
    # minimal path length along the mycelium L for anastomosis as in MycorrhizalPlant::simulateAnastomosis: anaMinDist if >= 0,
    # else distTH / sin(theta) + 2 dx (of the random parameter of the sub type)
    rp = plant.getOrganRandomParameter(ot, int(subtype))
    return rp.anaMinDist if rp.anaMinDist >= 0 else rp.distTH / max(np.sin(rp.theta), 1.e-6) + 2 * rp.dx

def getTipDistances(times, plant, rings, R_max = 0.1, ot = 5, exclude = "path"):
    # smallest distance [cm] of every actively growing tip to a hypha it could fuse with, at times[1:] (columns as getParasPerRing),
    # reconstructed from the final plant: the tips are those of getActiveTipsPerRingPerSubtype, and only segments created until t count.
    # The distance is the one of the anastomosis check of the model (MycorrhizalPlant::simulateAnastomosis, SDF_RootSystem::getDist):
    # distance to the segment minus its radius. inf if there is no such segment within R_max.
    # exclude = "path" (as the model): segments closer than L along the mycelium are excluded (anastomosisMinPath), fusion within
    #           the own hyphal tree is possible further away; "tree" (older version): all segments of the own hyphal tree are excluded.
    # With distTH = R a tip with dist < R would fuse at this check, so the fraction of tips with dist < R can be evaluated afterwards for any R <= R_max
    # (the model checks once per time step, against the segments of the step before, this uses the segments until t, the difference is small).
    # returns a dict with one entry per tip in a ring and time: "col" (index into times[1:]), "ring", "subtype" and "dist" [cm]
    t_eval = np.asarray(times)[1:]
    t_prev = np.asarray(times)[:-1]
    # segments of the hyphae (filter without pack keeps the node indices)
    ana = pb.SegmentAnalyser(plant)
    ana.filter("organType", ot)
    nodes = np.array([[n.x, n.y, n.z] for n in ana.nodes])
    segs = np.array([[s.x, s.y] for s in ana.segments], dtype=int).reshape(-1, 2)
    seg_ct = np.array(ana.getParameter("creationTime"))
    seg_tree = np.array(ana.getParameter("hyphalTreeIndex")).astype(int)
    seg_radius = np.array(ana.getParameter("radius"))
    a, b = nodes[segs[:, 0]], nodes[segs[:, 1]]
    mid = 0.5 * (a + b)
    reach = R_max + (0.5 * np.linalg.norm(b - a, axis=1).max() + seg_radius.max() if len(segs) else 0.) # segment midpoints that can be within R_max
    # growing tips at each time (as in getActiveTipsPerRingPerSubtype)
    col, ring_, subtype, tree, pos = [], [], [], [], []
    tip_organ, tip_k = [], []                                       # organ and local node index of each tip (exclude = "path")
    ringOfNode = {}
    for o in plant.getOrgans(ot):
        n = o.getNumberOfNodes()
        cts = np.array([o.getNodeCT(i) for i in range(n)])
        active = (t_eval >= cts[0]) & (t_prev < cts[-1])
        if not active.any():
            continue
        st, tr = int(o.getParameter("subType")), int(o.getParameter("hyphalTreeIndex"))
        tipIdx = np.searchsorted(cts, t_eval[active], side="right") - 1
        for j, k in zip(np.nonzero(active)[0], tipIdx):
            nodeId = o.getNodeId(int(k))
            if nodeId not in ringOfNode:
                ringOfNode[nodeId] = -1
                p = o.getNode(int(k))
                for r, ring in enumerate(rings):
                    if ring.getDist(p) < 0:
                        ringOfNode[nodeId] = r
                        break
            if ringOfNode[nodeId] >= 0:
                p = o.getNode(int(k))
                col.append(j), ring_.append(ringOfNode[nodeId]), subtype.append(st), tree.append(tr), pos.append([p.x, p.y, p.z])
                tip_organ.append(o), tip_k.append(int(k))
    col, tree, pos = np.array(col, dtype=int), np.array(tree, dtype=int), np.array(pos).reshape(-1, 3)
    dist = np.full(len(col), np.inf)
    if exclude == "path":
        organ_by_id = {o.getId(): o for o in plant.getOrgans(ot)}
        seg_oid = np.array(ana.getParameter("id")).astype(int)     # organ of each segment
        L_st = {}                                                   # minimal path length per sub type
        chains = {}                                                 # (organ id, local node index, L) -> chain
        local = {}                                                  # organ id -> {global node id: local index}
        def chain_of(o, k, L):
            key = (o.getId(), int(k), L)
            if key not in chains:
                chains[key] = _chain(o, k, L)
            return chains[key]
        def local_index(o, node):
            oid = o.getId()
            if oid not in local:
                local[oid] = {o.getNodeId(i): i for i in range(o.getNumberOfNodes())}
            return local[oid].get(node, -1)
    for j in np.unique(col):
        tips = np.nonzero(col == j)[0]
        existing = np.nonzero(seg_ct <= t_eval[j] + 1.e-9)[0]
        if len(existing) == 0:
            continue
        cand = cKDTree(mid[existing]).query_ball_point(pos[tips], reach)
        lens = np.array([len(c) for c in cand])
        if lens.sum() == 0:
            continue
        ti = np.repeat(tips, lens)
        si = existing[np.concatenate([c for c in cand if c]).astype(int)]
        if exclude == "tree":
            other = seg_tree[si] != tree[ti]
            ti, si = ti[other], si[other]
            d = _pointSegmentDistance(pos[ti], a[si], b[si]) - seg_radius[si]
            np.minimum.at(dist, ti, d)
            continue
        # as the model: the nearest segment that is at least L away along the mycelium; candidates in order of distance,
        # usually only the first few (the own hypha) are excluded
        d = _pointSegmentDistance(pos[ti], a[si], b[si]) - seg_radius[si]
        order = np.lexsort((d, ti))
        ti, si, d = ti[order], si[order], d[order]
        start = np.searchsorted(ti, tips)
        stop = np.searchsorted(ti, tips, side="right")
        for t, i0, i1 in zip(tips, start, stop):
            o_t, k_t = tip_organ[t], tip_k[t]
            st = int(subtype[t])
            if st not in L_st:
                L_st[st] = anastomosisMinPath(plant, st, ot)
            L = L_st[st]
            cA = chain_of(o_t, k_t, L)
            for i in range(i0, i1):
                if d[i] > R_max:
                    break
                o_s = organ_by_id.get(int(seg_oid[si[i]]))
                if o_s is None:
                    continue
                k_s = local_index(o_s, segs[si[i], 1])
                if k_s < 0 or (o_s.getId() == o_t.getId() and k_s == k_t):
                    continue                                        # tip itself
                if _pathDist(cA, chain_of(o_s, k_s, L)) >= L:
                    dist[t] = d[i]
                    break
    dist[dist > R_max] = np.inf
    return {"col": col, "ring": np.array(ring_, dtype=int), "subtype": np.array(subtype, dtype=int), "dist": dist}

def tipDistanceCounts(sim, R_values, subtype = None):
    # from the tip distances of one simulation ("tipdist_*" of petridishsimauto.py): number of growing tips (rings x columns, as "active_tips")
    # and number of them closer than R to a hypha they could fuse with (see getTipDistances) (len(R_values) x rings x columns), optionally only tips of one sub type
    shape = np.shape(sim["length_densities"])
    col, ring, d = sim["tipdist_col"], sim["tipdist_ring"], sim["tipdist_d"]
    if subtype is not None:
        m = sim["tipdist_subtype"] == subtype
        col, ring, d = col[m], ring[m], d[m]
    n_tips = np.zeros(shape)
    np.add.at(n_tips, (ring, col), 1)
    n_within = np.zeros((len(R_values),) + shape)
    for i, R in enumerate(R_values):
        m = d < R
        np.add.at(n_within[i], (ring[m], col[m]), 1)
    return n_tips, n_within

def proportionalFit(rho, events, tips):
    # C of f = C rho with f = events / tips: weighted least squares through the origin, minimises sum tips (events/tips - C rho)^2
    den = np.sum(tips * rho ** 2)
    return np.sum(events * rho) / den if den > 0 else np.nan

def tipDistanceC(simulations, R_values, min_tips = 5, subtype = None):
    # C(R) of f = C rho for every R, f = fraction of the growing tips closer than R to a hypha they could fuse with, rho = filament density
    # of the ring [um mm^-2], over all rings, times (first column left out, see getParasPerRing) and simulations with at least min_tips tips.
    # Returns C (len(R_values)) and the pooled data "rho", "tips", "within" (len(R_values) x points)
    rho_all, tips_all, within_all = [], [], []
    for sim in simulations:
        n_tips, n_within = tipDistanceCounts(sim, R_values, subtype)
        rho = sim["length_densities"][:, 1:] * 1000
        ok = n_tips[:, 1:] >= min_tips
        rho_all.append(rho[ok]), tips_all.append(n_tips[:, 1:][ok]), within_all.append(n_within[:, :, 1:][:, ok])
    rho, tips, within = np.concatenate(rho_all), np.concatenate(tips_all), np.concatenate(within_all, axis=1)
    C = np.array([proportionalFit(rho, w, tips) for w in within])
    return C, {"rho": rho, "tips": tips, "within": within}

def encounterData(simulations):
    # realized anastomoses of the model per ring and interval: {"total", "runner", "bas": {"rho" [um mm^-2], "ana" (new anastomoses per hour),
    # "tips" (growing tips)}}, one entry per ring, interval and simulation, first column left out (see getParasPerRing)
    # (needs "ana_new*", "active_tips*" of petridishsimauto.py), f = ana / tips is the fraction of the growing tips that fuse per hour
    suffix = {"total": "", "runner": "_runner", "bas": "_bas"}
    data = {key: {"rho": [], "ana": [], "tips": []} for key in suffix}
    for sim in simulations:
        dt_h = np.diff(np.concatenate([[0.], sim["times"]]))[1:] * 24 # interval lengths [h]
        rho = sim["length_densities"][:, 1:] * 1000
        for key, s in suffix.items():
            ana, tips = sim["ana_new" + s][:, 1:], sim["active_tips" + s][:, 1:]
            data[key]["rho"].append(rho.ravel())
            data[key]["ana"].append((ana / dt_h).ravel())
            data[key]["tips"].append(tips.ravel())
    return {key: {k: np.concatenate(v) for k, v in d.items()} for key, d in data.items()}

def encounterC(simulations, min_tips = 5):
    # C of f = C rho for the realized anastomoses per hour (encounterData), {"total", "runner", "bas": C}, only points with at least min_tips tips
    result = {}
    for key, d in encounterData(simulations).items():
        ok = d["tips"] >= min_tips
        result[key] = proportionalFit(d["rho"][ok], d["ana"][ok], d["tips"][ok]) if ok.any() else np.nan
    return result

def anastomosisAlpha(simulations, min_tips = 5):
    # alpha [um h^-1] of the anastomosis rate a = alpha rho n. Nomenclature here: alpha = anastomosis, beta = branching
    # (branchingRate); Oyarte Galvez et al. 2025 use the two letters the other way round (their beta is this alpha, about 23 um/h).
    # a = new anastomoses per area and hour [mm^-2 h^-1],
    # rho = filament density of all hyphae [mm mm^-2], n = growing tip density [mm^-2], per ring and interval (first column left out,
    # see getParasPerRing), only points with at least min_tips growing tips. alpha = sum a / sum rho n over all points and simulations,
    # i.e. the scale between the two curves a(t) and rho n(t) of the paper. Per sub type: a and n of the tips of that sub type, rho of all
    # hyphae. The fraction of the growing tips that fuse per hour is f = a / n = alpha rho (alpha 23 um/h, rho 0.87 mm mm^-2: f = 2 %/h).
    # Same as C of encounterC (alpha = 1e6 C) up to the weighting of the points. Returns {"total", "runner", "bas": alpha}, nan without points
    suffix = {"total": "", "runner": "_runner", "bas": "_bas"}
    a_sum = dict.fromkeys(suffix, 0.)
    rn_sum = dict.fromkeys(suffix, 0.)
    for sim in simulations:
        dt_h = np.diff(np.concatenate([[0.], sim["times"]]))[1:] * 24 # interval lengths [h]
        area = np.asarray(sim["area"]).reshape(-1, 1)                 # ring areas [mm^2]
        rho = sim["length_densities"][:, 1:]                           # [mm mm^-2]
        for key, s in suffix.items():
            ana, tips = sim["ana_new" + s][:, 1:], sim["active_tips" + s][:, 1:]
            ok = tips >= min_tips
            a_sum[key] += np.sum((ana / dt_h / area)[ok])
            rn_sum[key] += np.sum((rho * tips / area)[ok])
    return {key: 1000 * a_sum[key] / rn_sum[key] if rn_sum[key] > 0 else np.nan for key in suffix}

anastomosisBeta = anastomosisAlpha # old name (letters as in Oyarte Galvez et al. 2025), kept for older scripts

def theoryC(R, layer, dim = 3):
    # C [mm^2 um^-1] of a random network of straight hyphae (Poisson line process, any orientation): the expected number of hyphae closer
    # than R to a point is pi R^2 rho_V in 3D (rho_V = rho / layer, length per volume) and 2 R rho in 2D (rho length per area),
    # for small R the fraction is f = C rho. R and layer [cm], rho [um mm^-2]
    R_mm, layer_mm = np.asarray(R) * 10, layer * 10
    return np.pi * R_mm ** 2 / layer_mm / 1000 if dim == 3 else 2 * R_mm / 1000

def sigmoidLengthDens(t, K1,lamda_dens,t_arrival):
    rho = K1/(1+np.exp(np.minimum(lamda_dens*(t_arrival-t), 700.))) # limited to avoid overflow, exp(700) is already infinite for the fit
    return rho

def sigmoidTipDens(t, K2,lamda_tips,t_arrival):
    # derivative of the sigmoid (bell curve), scaled such that K2 is the peak height at t = t_arrival
    e = np.exp(-np.abs(lamda_tips*(t_arrival-t))) # symmetric, using the smaller exponent avoids overflow
    n = 4 * K2 * e / (1+e)**2
    return n

def arrivalTimeDensFit(t,rho):
    # fits sigmoidLengthDens to the length density rho(t) of one ring, returns (K1, lamda_dens, t_arrival) and their covariance
    span = t[-1] - t[0]
    p0 = [max(rho.max(), 1.e-12), 4 / span, t[np.argmax(rho >= rho.max() / 2)]]
    bounds = ([0, 0, t[0] - span], [np.inf, np.inf, t[-1] + span])
    param, param_cov = curve_fit(sigmoidLengthDens, t, rho, p0, bounds=bounds, maxfev=10000)
    return param, param_cov

def arrivalTimeTipFit(t, n, t_arrival, lamda):
    # fits sigmoidTipDens to the tip density n(t) of one ring, starting from the arrival time and rate of the length density fit
    span = t[-1] - t[0]
    p0 = [max(n.max(), 1.e-12), lamda, t_arrival]
    bounds = ([0, 0, t[0] - span], [np.inf, np.inf, t[-1] + span])
    param, param_cov = curve_fit(sigmoidTipDens, t, n, p0, bounds=bounds, maxfev=10000)
    return param, param_cov

def fitArrivalTimes(t, rho, tips, min_fraction=0.5):
    # arrival time of the network at each ring (rows of rho and tips, columns at times t) from a sigmoid fit of the length density rho,
    # tips are fitted with the sigmoid derivative. Only rings where the arrival lies within t and rho reaches min_fraction of the fitted
    # saturation K1 are kept. Returns a list of dicts with "ring", "t_n", "length" (K1, lamda, t_arrival) and "tips" (K2, lamda, t_arrival or None),
    # and the coefficient of determination of both fits "R2_length", "R2_tips" (nan without tip fit)
    def R2(y, y_fit):
        ss = np.sum((y - y.mean())**2)
        return 1 - np.sum((y - y_fit)**2) / ss if ss > 0 else np.nan
    fits = []
    for k in range(rho.shape[0]):
        if rho[k].max() <= 0:
            continue
        try:
            pl, _ = arrivalTimeDensFit(t, rho[k])
        except (RuntimeError, ValueError):
            continue
        if not (t[0] < pl[2] < t[-1]) or rho[k, -1] < min_fraction * pl[0]:
            continue
        try:
            pt, _ = arrivalTimeTipFit(t, tips[k], pl[2], pl[1])
        except (RuntimeError, ValueError):
            pt = None
        fits.append({"ring": k, "t_n": pl[2], "length": pl, "tips": pt, "R2_length": R2(rho[k], sigmoidLengthDens(t, *pl)),
                     "R2_tips": R2(tips[k], sigmoidTipDens(t, *pt)) if pt is not None else np.nan})
    return fits

FWHM_BELL = 2 * np.log(3 + 2 * np.sqrt(2)) # full width at half maximum of sigmoidTipDens is FWHM_BELL / lamda (about 3.53 / lamda)

def fitSpatialFronts(simulations, location, min_R2 = 0.5, min_rings = 5):
    # the tip wave in space: at every time of every simulation the tip density over the rings n(r) is fitted with the bell curve
    # sigmoidTipDens(r, A, k, r_f) (same function as over time): r_f = position of the front (peak of the tips) [mm], A = peak tip density
    # (amplitude) [mm^-2], width = FWHM_BELL / k = full width at half maximum of the tip band [mm]. Kept only if the peak and both
    # half maximum points r_f +- width/2 lie within the ring centers (else amplitude and width are extrapolated), R2 >= min_R2
    # and at least min_rings rings have tips.
    # simulations as from load_petridish_results ("times" days since t*, "tip_dens" [mm^-2]), location: ring centers [cm];
    # first column left out (see getParasPerRing). Returns a list of dicts "seed", "j" (column of tip_dens[:, 1:]), "t" [h], "r_f", "A", "k", "width", "R2"
    r = np.asarray(location) * 10
    fits = []
    for sim in simulations:
        t_h = sim["times"][1:] * 24
        n = sim["tip_dens"][:, 1:]
        for j, t in enumerate(t_h):
            y = n[:, j]
            if np.count_nonzero(y) < min_rings or y.max() <= 0:
                continue
            span = r[-1] - r[0]
            p0 = [y.max(), 8 / span, r[np.argmax(y)]]
            bounds = ([0, 0, r[0] - span], [np.inf, np.inf, r[-1] + span])
            try:
                p, _ = curve_fit(sigmoidTipDens, r, y, p0, bounds=bounds, maxfev=10000)
            except (RuntimeError, ValueError):
                continue
            ss = np.sum((y - y.mean())**2)
            R2 = 1 - np.sum((y - sigmoidTipDens(r, *p))**2) / ss if ss > 0 else np.nan
            half = FWHM_BELL / p[1] / 2
            if not (r[0] <= p[2] - half and p[2] + half <= r[-1]) or not (R2 >= min_R2): # peak and both half maximum points within the rings
                continue
            fits.append({"seed": int(sim["seed"]), "j": j, "t": t, "r_f": p[2], "A": p[0], "k": p[1], "width": FWHM_BELL / p[1], "R2": R2})
    return fits

def theilSlope(x, y):
    # robust slope and intercept of y against x (Theil-Sen, median of the pairwise slopes), nan for less than 3 points
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if len(x) < 3 or np.unique(x).size < 2:
        return np.nan, np.nan
    from scipy.stats import theilslopes
    return tuple(theilslopes(y, x)[:2])

def seedSummary(name, x, unit = "", test_zero = False, n_boot = 10000, seed = 1, verbose = True):
    # one value per seed (the seeds are the independent replicates): mean, median, sd, coefficient of variation and 95% bootstrap
    # confidence interval of the mean over the seeds, optionally a t-test against 0. Returns a dict (None for less than 2 seeds)
    from scipy.stats import ttest_1samp
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        if verbose:
            print(f"  {name}: less than 2 seeds")
        return None
    rng = np.random.default_rng(seed)
    lo, hi = np.percentile(x[rng.integers(0, len(x), (n_boot, len(x)))].mean(axis=1), [2.5, 97.5])
    res = {"mean": x.mean(), "median": np.median(x), "sd": x.std(ddof=1), "cv": x.std(ddof=1) / abs(x.mean()) if x.mean() != 0 else np.nan,
           "lo": lo, "hi": hi, "p": ttest_1samp(x, 0).pvalue if test_zero else np.nan, "n": len(x)}
    if verbose:
        print(f"  {name}: mean {res['mean']:.4g} {unit}, median {res['median']:.4g}, sd {res['sd']:.3g} (cv {res['cv']:.2f}), "
              f"95% CI of mean [{lo:.4g}, {hi:.4g}], n = {len(x)} seeds" + (f", t-test mean = 0: p = {res['p']:.3g}" if test_zero else ""))
    return res

def collectArrivalFits(simulations):
    # fitArrivalTimes for all simulations (dicts with "times" [days since t*], "length_densities" [mm mm^-2], "tip_dens" [mm^-2]),
    # time in h and length density in um mm^-2. Adds "seed" and the data "t", "rho_data", "tip_data" to each fit
    # the first column of length_densities only counts the first interval (see getParasPerRing), so it is left out
    fits = []
    for sim in simulations:
        t_h = sim["times"][1:] * 24
        rho = sim["length_densities"][:, 1:] * 1000
        tips = sim["tip_dens"][:, 1:]
        for f in fitArrivalTimes(t_h, rho, tips):
            f.update({"seed": int(sim["seed"]), "t": t_h, "rho_data": rho[f["ring"]], "tip_data": tips[f["ring"]]})
            fits.append(f)
    return fits
