import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import plantbox as pb
from scipy.optimize import curve_fit

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
