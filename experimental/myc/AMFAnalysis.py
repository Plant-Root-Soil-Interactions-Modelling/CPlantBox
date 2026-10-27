import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import plantbox as pb
import time
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import plantbox.visualisation.vtk_plot as vp
import networkx as nx

def getMycSegmentAnalyser(plant, write=[],filename = "write",step = 0, stdwrite = False,):
        ana = pb.SegmentAnalyser(plant)
        ana.addData("colonization", plant.getNodeColonizations(2))
        ana.addData("colonizationTime", plant.getNodeColonizationTime(2))
        ana.addData("anastomosis", plant.getAnastomosisPoints(5))
        ana.addData("nodeTips", plant.getNodeTips(5))
        # ana.crop(crop)
        if write:
            ana.write(filename + "_" + str(step) + "_" +time.strftime("%Y%m%d-%H%M%S") + ".vtp",write)
        if stdwrite:
            ana.write(filename + "_" + str(step) + "_" +time.strftime("%Y%m%d-%H%M%S") + ".vtp",["radius","creationTime","subType","organType","colonization","colonizationTime","anastomosis","nodeTips"])
        return ana


def getLengthPerSubtype(plant):
    hyphae = pb.SegmentAnalyser(plant)
    hyphae.filter("organType",5.0)
    hyphae.pack()
    lenSubtype = []
    for i in range(1,4):
        hyphae.filter("subType",i)
        hyphae.pack()
        lenSubtype.append(hyphae.getSummed("length"))
        hyphae = pb.SegmentAnalyser(plant)
        hyphae.filter("organType",5.0)
        hyphae.pack()
    return lenSubtype

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

def getParasPerRing(parameters, times, plant, rings):
    # same matrices as getParaSumperRing (rings x len(times)-1), but for all parameters at once:
    # one SegmentAnalyser, ring of each segment computed once, time intervals binned with numpy
    ana = getMycSegmentAnalyser(plant)
    ringIdx = getRingIndex(ana, rings)
    ct = np.array(ana.getParameter("creationTime"))
    n = len(times) - 1
    timeIdx = np.searchsorted(times, ct, side="left") # 0: ct <= times[0], i: times[i-1] < ct <= times[i]
    valid = (ringIdx >= 0) & (ct >= 0) & (timeIdx <= n)
    mats = {}
    for parameter in parameters:
        values = np.array(ana.getParameter(parameter))
        summed = np.zeros((len(rings), n + 1))
        np.add.at(summed, (ringIdx[valid], timeIdx[valid]), values[valid])
        if parameter == "nodeTips": # we want to look only at actively growing tips. (Oyarte Galvez 2025 page 3 when discussing Fig 2b)
            mat = summed[:, 1:] # created in (times[i], times[i+1]]
        else:
            mat = np.cumsum(summed, axis=1)[:, 1:] # created until times[i+1]
            mat[:, 0] = summed[:, 1] # first column only counts (times[0], times[1]], as in getParaSumperRing
        mats[parameter] = mat
    return mats

def getParaSumperRing(parameter, times, plant, rings):
    return getParasPerRing([parameter], times, plant, rings)[parameter]

def getActiveTipsPerRing(times, plant, rings, ot = 5):
    # number of actively growing tips per ring at times[1:] (rings x len(times)-1, same layout as getParasPerRing),
    # reconstructed from the final plant: an organ is growing at time t if it was created (node 0) before t and
    # its last node after t (or it is still active at the end). Its tip at time t is its last node created until t.
    # (node 1 is not used as start, it gets a new creation time as long as the first segment is shorter than dx and the node is moved)
    # ("nodeTips" in getParasPerRing instead counts organs whose last node was created in the interval, i.e. that stopped growing then)
    t_eval = np.asarray(times)[1:]
    mat = np.zeros((len(rings), len(t_eval)))
    ringOfNode = {} # global node id -> ring index (-1 = in no ring)
    for o in plant.getOrgans(ot):
        n = o.getNumberOfNodes()
        cts = np.array([o.getNodeCT(i) for i in range(n)])
        t_start, t_stop = cts[0], cts[-1] # node 0 is the branching point on the parent, its creation time is the emergence of the organ
        active = (t_eval >= t_start) & ((t_eval < t_stop) | o.isActive())
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
    return mat

def getParaSumperRing_old(parameter, times, plant, rings):
    # print("Parameter sum for "+ parameter)
    paradenmat = np.zeros((len(rings),len(times[1:])))
    flipped = np.flip(np.asarray(times))
    for k, ring in enumerate(rings):
        ana = getMycSegmentAnalyser(plant)
        for j in range(len(times[1:])-1):
            ana.filter("creationTime", 0, flipped[j])
            ana.pack()
            distrib = ana.getSummed(parameter, ring)
            ana.filter("creationTime",0,flipped[j+1])
            ana.pack()
            summed = 0
            if parameter == "nodeTips": # we want to look only at actively growing tips. (Oyarte Galvez 2025 page 3 when discussing Fig 2b)
                summed = ana.getSummed(parameter, ring)
            paradenmat[k, len(times[1:])-1 -j] = np.array(distrib-summed).sum()
        ana.filter("creationTime",0,flipped[len(times[1:])-1])
        ana.pack()
        summed = ana.getSummed(parameter, ring)
        ana.filter("creationTime",0,flipped[len(times[1:])])
        ana.pack()
        summed = summed - ana.getSummed(parameter, ring)
        paradenmat[k, 0] = np.array(summed).sum()
    return paradenmat

def getParameterOverTime(parameter, times, plant, subType):
    paradenmat = np.zeros((len(subType), len(times[1:])))
    flipped = np.flip(np.asarray(times))
    for k in subType:
        ana = pb.SegmentAnalyser(plant)
        ana.filter("organType", 5, 5)
        for j in range(len(times[1:])-1):
            ana.filter("creationTime", 0, flipped[j])
            ana.pack()
            distrib = ana.getSummed(parameter)
            ana.filter("creationTime",0,flipped[j+1])
            ana.pack()
            summed = ana.getSummed(parameter)
            paradenmat[k-1,len(times[1:])-1 -j] = np.array(distrib-summed).sum() #/np.array(allTips-overTips).sum() if np.array(allTips-overTips).sum() != 0 else 0
        ana.filter("creationTime",0,flipped[len(times[1:])-1])
        ana.pack()
        summed = ana.getSummed(parameter)
        ana.filter("creationTime",0,flipped[len(times[1:])])
        ana.pack()
        summed = summed - ana.getSummed(parameter)
        paradenmat[k-1,-1] = np.array(summed).sum()
    return paradenmat

def makedishes(diameter, height, barrier_thickness, barrier_height, opening_length, opening_height, rootradius,nRings):
    radius = diameter / 2
    # petri dish has a radius of 9.4 cm and a height of 1.6 cm
    petri_dish = pb.SDF_PlantContainer(radius,radius,height,False)
    # the helper dish is used to cut the petri dish in half, it has the same radius and height as the petri dish but is rotated and translated to cut the petri dish in half
    helper_dish = pb.SDF_PlantContainer(radius,radius,height,True)
    # moving the helper dish such that it halves the petri dish and removes a bit more to restrict roots to the correct side of barrier
    moved_helper_dish = pb.SDF_RotateTranslate(helper_dish, 0, 0, pb.Vector3d(-(radius+barrier_thickness+rootradius), 0, 0))
    half_dish = pb.SDF_Intersection(petri_dish, moved_helper_dish)
    # helper container for barrier
    helper_staff = pb.SDF_PlantBox(barrier_thickness,barrier_height,diameter)
    # helper container for opening in barrier
    helper_staff2 = pb.SDF_PlantBox(barrier_thickness,opening_height,opening_length)
    # since master (3b88115e) SDF_PlantBox is shifted down by eps = 0.25 cm, which would move barrier and opening; keep the old geometry
    if hasattr(helper_staff, "eps"):
        helper_staff.eps = 0.
        helper_staff2.eps = 0.
    # have to  move the helper staff for the right position of the opening and barrier, the midpoint is the distance from the center of the helper staff to the center of the petri dish, the bottompoint is the distance from the center of the helper staff to the center of the petri dish in the y direction, and the helper staff is moved to the position of the opening and barrier
    midpoint = opening_length / 2 - diameter / 2
    bottompoint = opening_height / 2 - barrier_height / 2
    # container for opening moved to the right position
    helper_staff2 = pb.SDF_RotateTranslate(helper_staff2, 0,0,pb.Vector3d(0, bottompoint, midpoint))
    # opening made in the barrier
    helper_dish2 = pb.SDF_Difference(helper_staff, helper_staff2)
    # barrier moved to the right position
    moved_helper_dish2 = pb.SDF_RotateTranslate(helper_dish2, 90, pb.SDF_Axis.xaxis , pb.Vector3d(0, -radius, -height/2))
    petri_dish = pb.SDF_Difference(petri_dish, moved_helper_dish2)

    moved_helper_dish_hyphae = pb.SDF_RotateTranslate(helper_dish, 0, 0, pb.Vector3d(-(radius+barrier_thickness/2), 0, 0))
    small_dish = pb.SDF_PlantContainer(radius,radius,height,False)
    small_hyphae_dish = pb.SDF_Difference(small_dish, moved_helper_dish_hyphae)

    small_dish = pb.SDF_PlantContainer(radius*np.sqrt(1/nRings),radius*np.sqrt(1/nRings),height,False)
    ringone = pb.SDF_Difference(small_dish, moved_helper_dish_hyphae)
    rings = []
    rings.append(ringone)
    for i in range(2, nRings+1):
        small_dish = pb.SDF_PlantContainer(radius*np.sqrt(i/nRings),radius*np.sqrt(i/nRings),height,False)
        small_hyphae_dish2 = pb.SDF_Difference(small_dish, moved_helper_dish_hyphae)
        old_dish = pb.SDF_Difference(pb.SDF_PlantContainer(radius*np.sqrt((i-1)/nRings),radius*np.sqrt((i-1) /nRings),height,False),moved_helper_dish_hyphae)
        small_hyphae_dish2 = pb.SDF_Difference(small_hyphae_dish2,old_dish)
        rings.append(small_hyphae_dish2)
    return petri_dish, small_hyphae_dish, half_dish, rings

def makeBoxes(radius, height, boxlength, barrier_thickness, small_hyphae_dish, resolution = 0.005):
    # non overlapping square boxes (edge length boxlength) tiling the hyphae side of the dish (x > -barrier_thickness/2),
    # each box is cut with small_hyphae_dish, boxes outside the dish are dropped
    # returns boxes, centres of the uncut boxes (n x 2) and volume of each cut box
    box = pb.SDF_PlantContainer(boxlength/2, boxlength/2, height, True) # square container has half edge length r
    x_left = -barrier_thickness/2
    nx = int(np.ceil((radius - x_left)/boxlength))
    ny = int(np.ceil(2*radius/boxlength))
    # area of the cut boxes on a fine grid, same geometry as small_hyphae_dish
    s = (np.arange(round(boxlength/resolution)) + 0.5) * resolution - boxlength/2
    sx, sy = np.meshgrid(s, s)
    boxes, centres, volumes = [], [], []
    for i in range(nx):
        for j in range(ny):
            cx = x_left + (i + 0.5)*boxlength
            cy = -ny*boxlength/2 + (j + 0.5)*boxlength
            x, y = sx + cx, sy + cy
            area = np.count_nonzero((x*x + y*y < radius*radius) & (x > x_left)) * resolution*resolution
            if area > 0:
                moved_box = pb.SDF_RotateTranslate(box, 0, 0, pb.Vector3d(cx, cy, 0))
                boxes.append(pb.SDF_Intersection(moved_box, small_hyphae_dish))
                centres.append([cx, cy])
                volumes.append(area*height)
    return boxes, np.array(centres), np.array(volumes)

def ringLocations(radius, nRings):
    # outer radius of each ring, rings have equal area (see makedishes)
    return np.array([radius*np.sqrt(i/nRings) for i in range(1, nRings+1)])

def ringCenters(outer):
    # radius that splits each ring into two halves of equal area, from the outer radii (ringLocations)
    inner = np.concatenate(([0.], outer[:-1]))
    return np.sqrt((inner*inner + outer*outer) / 2)

def ringVolume(radius, nRings, height, barrier_thickness = 0.):
    # volume of each ring (array of length nRings), the rings are cut at x = -barrier_thickness/2 (see makedishes),
    # so besides the half annulus they contain the strip -barrier_thickness/2 < x < 0
    c = barrier_thickness / 2
    R = np.concatenate(([0.], ringLocations(radius, nRings)))
    c_ = np.minimum(c, R)
    # area of the disk of radius R with x > -c: half disk + segment between x = -c and x = 0
    area = np.pi / 2 * R*R + c_*np.sqrt(R*R - c_*c_) + R*R*np.arcsin(np.divide(c_, R, out=np.zeros_like(R), where=R > 0))
    return np.diff(area) * height

def setUpSimulationTime(simTime, fps):
    # set up simulation time
    simTime = simTime
    fps = fps
    dt = 1 / fps
    nSteps = int(simTime / dt)
    return dt, nSteps

def makeSimulation(seed, path, name, height, petri_dish, small_hyphae_dish, half_dish, dt, nSteps, hours_hyphae, filename, animation = False, verbose = False,
                   max_wait_days = None, stats = None):
    # set up simulation
    # max_wait_days: stop waiting for t* (10 cm hyphae behind the barrier) after this many simulated days, then returns mycp, None, None
    #                (None = wait without limit, as before). A limit in simulated time is reproducible, a wall clock limit is not.
    # stats: optional dict, filled with "status" ("ok" or "not_crossed"), the wall clock time of each phase [s] ("wall_root", "wall_wait",
    #        "wall_hyphae"), "t_root_end" (simulation time after the root phase), "wait_days" (simulated days until t*, or until the stop)
    #        and the length behind the barrier after each waiting step ("crossed_t", "crossed_len")
    #### WORK IN PROGRESS
    t_wall = time.perf_counter()
    mycp = pb.MycorrhizalPlant(seed)
    mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)
    seed_para = pb.SeedRandomParameter(mycp)
    seed_para.seedPos.z = -height / 2
    seed_para.seedPos.x = -9.4 / 2 + 0.5
    seed_para.seedPos.y = 0
    mycp.setOrganRandomParameter(seed_para)
    mycp.setGeometry(half_dish)
    mycp.initialize()

    # print("Starting simulation with seed: " + str(seed))
    for i in range(nSteps):
        # if (i % 10 == 0):
        #     print("Step " + str(i) + " of " + str(nSteps))
        mycp.simulate(dt,verbose)
        if (animation):
            getMycSegmentAnalyser(mycp,filename = "animation/" + filename + "_anim",step = i, stdwrite= True)    # look at roots and container

    root = mycp.getOrganRandomParameter(pb.root)
    for rp in root:
        rp.hyphalEmergenceDensity = 4.0
        rp.lmbd = 0.15
        mycp.setOrganRandomParameter(rp)

    mycp.changeGeometry(5,petri_dish)

    wall_root = time.perf_counter() - t_wall
    t_wall = time.perf_counter()
    t_root_end = mycp.getSimTime()
    crossed_t, crossed_len = [], []
    status = "ok"
    crossed_barrier = 0

    while crossed_barrier < 10:
        if max_wait_days is not None and mycp.getSimTime() - t_root_end >= max_wait_days - 1.e-9:
            status = "not_crossed"
            break
        mycp.simulate(dt,verbose)
        ana = pb.SegmentAnalyser(mycp) # only the length is needed, cheaper than getMycSegmentAnalyser (no colonization, anastomosis, tips)
        ana.filter("organType",5)
        ana.pack()
        crossed_barrier = ana.getSummed("length",small_hyphae_dish)
        crossed_t.append(mycp.getSimTime())
        crossed_len.append(crossed_barrier)

    if stats is not None:
        stats.update({"status": status, "wall_root": wall_root, "wall_wait": time.perf_counter() - t_wall, "wall_hyphae": 0.,
                      "t_root_end": t_root_end, "wait_days": mycp.getSimTime() - t_root_end,
                      "crossed_t": np.array(crossed_t), "crossed_len": np.array(crossed_len)})
    if status != "ok":
        print("Seed " + str(seed) + ": less than 10 cm hyphae behind the barrier after " + str(max_wait_days) + " days, stopped")
        return mycp, None, None
    t_wall = time.perf_counter()

    time_star = mycp.getSimTime()
    for i in range(0, hours_hyphae):
        # if (i % 10 == 0):
        #     print("Step " + str(i) + " of " + str(hours_hyphae))
        mycp.simulate(dt,verbose)
        if (animation):
            getMycSegmentAnalyser(mycp,filename = "animation/" + filename + "_anim",step = i+nSteps, stdwrite= True)    # look at roots and container
    if stats is not None:
        stats["wall_hyphae"] = time.perf_counter() - t_wall
    nmbtimes = (max(mycp.getParameter("creationTime")) - time_star)*24
    nmbtimespoint = (max(mycp.getParameter("creationTime")) + 0.1 - time_star)*24
    times = np.linspace(time_star, max(mycp.getParameter("creationTime"))+0.01, max(round(nmbtimes),round(nmbtimespoint)))
    # times = times - time_star
    return mycp, time_star, times

def sigmoidLengthDens(t, K1,lamda_dens,t_arrival):
    rho = K1/(1+np.exp(lamda_dens*(t_arrival-t)))
    return rho

def sigmoidTipDens(t, K2,lamda_tips,t_arrival):
    n = (K2 * np.exp(lamda_tips*(t_arrival-t)))/(1+np.exp(lamda_tips*(t_arrival-t)))
    return n

def arrivalTimeDensFit(t,rho):
    p0 = [max(rho), np.median(t),1, min(rho)]
    param, param_cov = curve_fit(sigmoidLengthDens, t, rho, p0)
    return param, param_cov

def prepare_data(lenMat, vol, anaMat, tipMat):
    rho = lenMat / np.reshape(vol, (-1, 1)) # vol: one value or one per ring/box
    # print(rho)
    f_enc = np.divide(
        anaMat, tipMat,
        out=np.zeros_like(anaMat, dtype=float),
        where= (tipMat != 0)
    )
    n_locations, n_times = rho.shape
    rho = rho.ravel()
    f_enc = f_enc.ravel()
    # Zeile = Ort, Spalte = Zeit
    location = np.repeat(np.arange(n_locations), n_times)
    time = np.tile(np.arange(n_times), n_locations)
    valid = np.isfinite(rho) & np.isfinite(f_enc)
    return rho[valid], f_enc[valid], location[valid], time[valid]


def find_linear_range(rho, f_enc, window, step, min_points=10):
    results = []
    starts = np.arange(rho.min(), rho.max() - window, step)
    for start in starts:
        end = start + window
        mask = (rho >= start) & (rho <= end)
        if mask.sum() < min_points:
            continue
        r = rho[mask]
        f = f_enc[mask]
        model = LinearRegression(fit_intercept=False).fit(r.reshape(-1, 1), f)
        results.append({
            "rho_min": start,
            "rho_max": end,
            "rho_center": (start + end) / 2,
            "C": model.coef_[0],
            "R2": model.score(r.reshape(-1, 1), f),
            "n": mask.sum()
        })
    return results

def analyze_simulations(simulations, window, step, min_points=10):
    all_results = {}
    for sim in simulations:
        rho, f_enc, location, time = prepare_data(
            sim["lenMat"], sim["vol"], sim["anaMat"], sim["tipMat"]
        )
        all_results[sim["label"]] = {
            "rho": rho,
            "f_enc": f_enc,
            "location": location,
            "time": time,
            "linear": find_linear_range(
                rho, f_enc, window, step, min_points
            )
        }
    return all_results

def plot_data(all_results):
    fig, ax = plt.subplots(figsize=(8, 5.5))
    n_locations = max(
        r["location"].max() + 1
        for r in all_results.values()
    )
    cmap = plt.cm.viridis
    norm = plt.Normalize(0, n_locations - 1)
    markers = ["o", "s", "^", "D", "v", "P", "X"]
    for i, (label, data) in enumerate(all_results.items()):
        ax.scatter(
            data["rho"], data["f_enc"],
            c=data["location"], cmap=cmap, norm=norm,
            marker=markers[i % len(markers)],
            alpha=0.3, s=20, edgecolors="none",
            label=label
        )

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("location")

    ax.set_xlabel(r"realized hyphal length density $\rho_h$")
    ax.set_ylabel(r"$f_{\mathrm{enc}}$")
    ax.set_title(r"Encounter fraction vs. $\rho_h$")
    ax.grid(alpha=0.2)
    ax.legend()

    plt.tight_layout()
    plt.show()

def plot_linear_analysis(all_results):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    for label, data in all_results.items():
        linear = data["linear"]
        if not linear:
            continue
        rho = np.array([r["rho_center"] for r in linear])
        C = np.array([r["C"] for r in linear])
        R2 = np.array([r["R2"] for r in linear])
        ax1.plot(rho, C, "o-", ms=4, label=label)
        ax2.plot(rho, R2, "o-", ms=4, label=label)
    ax1.set_xlabel(r"density $\rho_h$")
    ax1.set_ylabel(r"$C$")
    ax1.set_title(r"Linear coefficient $C$")
    ax2.set_xlabel(r"density $\rho_h$")
    ax2.set_ylabel(r"$R^2$")
    ax2.set_title(r"Linearity")
    for ax in (ax1, ax2):
        ax.grid(alpha=0.2)
        ax.legend()
    plt.tight_layout()
    plt.show()

def interpolate_to_common_time(simulations, key, n_times=99):
    # interpolates sim[key] (rings x times) of all simulations onto the time range they share
    from scipy.interpolate import interp1d
    t_start = max(sim["times"][0] for sim in simulations)
    t_end = min(sim["times"][-1] for sim in simulations)
    if t_start >= t_end:
        raise ValueError(f"interpolate_to_common_time: time ranges of the simulations do not overlap (latest start {t_start}, earliest end {t_end})")
    t_common = np.linspace(t_start, t_end, n_times)
    interp = np.stack([interp1d(sim["times"], sim[key], axis=1)(t_common) for sim in simulations])
    return t_common, interp   # (n_sim, n_r, n_t)

def plot_ring_heatmap(times, location, mat, cbar_label, filename):
    plt.pcolormesh(times, location, mat, shading='auto', cmap='plasma')
    plt.colorbar(label=cbar_label)
    plt.xlabel("Time [days]")
    plt.ylabel("Distance [cm] from centre")
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()

def plot_radial_profiles(location, t_common, interp, indices, ylabel="Tip amount", stat="median", filename="TipDensOverAllSim.png"):
    # profiles over the simulations at the time indices given, location should be the ring centers (ringCenters)
    # stat = "median": median with 25%-75% quantiles (robust, the distribution over seeds is skewed), "mean": mean +- std
    if stat == "median":
        mid = np.median(interp, axis=0)
        low, high = np.percentile(interp, [25, 75], axis=0)
        band = "25%-75% quantile"
    else:
        mid = interp.mean(axis=0)
        std = interp.std(axis=0)
        low, high = mid - std, mid + std
        band = r"mean $\pm$ std"
    fig, ax = plt.subplots(figsize=(6, 5))
    colors = plt.get_cmap("plasma")(np.linspace(0.2, 0.9, len(indices)))
    for color, i in zip(colors, indices):
        ax.plot(location, mid[:, i], color=color, lw=2, marker="o", ms=3, label=f"t = {t_common[i]:.2f} d")
        ax.fill_between(location, low[:, i], high[:, i], color=color, alpha=0.3)
    ax.set_xlabel(r"Distance from centre $r$ (cm)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{stat} over {interp.shape[0]} simulations, band: {band}", fontsize=9)
    ax.legend(title=r"time after $t^*$")
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)

def save_simulation(filename, sim):
    # stores one simulation dict (only arrays/numbers, e.g. "times", "tip_dens", ...) as .npz
    np.savez(filename, **sim)

def load_simulations(pattern):
    # loads all .npz files matching pattern (e.g. "results/seed_*.npz") into a list of dicts, sorted by filename
    import glob
    simulations = []
    for f in sorted(glob.glob(pattern)):
        with np.load(f) as data:
            simulations.append({key: data[key] for key in data.files})
    # print("Loaded " + str(len(simulations)) + " simulations from " + pattern)
    return simulations


def estimate_carbon_cost_per_Org(length,radius,pdry,pcarbon,dcell,CUE):
    gamma = np.pi * radius *radius*dcell*pdry*pcarbon
    gamma = gamma / CUE
    carbon_cost = gamma*length
    return carbon_cost

def getAnastomosisMerges(mycp):
    # {tip node id: node id it merged with}, at an anastomosis the hypha gets a new node at the position
    # of the target node (see MycorrhizalPlant::simulateAnastomosis), both are the same point in the network
    merges = {}
    for h in mycp.getOrgans(pb.hyphae):
        mergeID = int(h.getParameter("mergePointID"))
        if mergeID > -1:
            merges[h.getNodeId(h.getNumberOfNodes() - 1)] = mergeID
    for tip in merges: # a target might itself be a merged tip
        target = merges[tip]
        while target in merges and merges[target] != tip:
            target = merges[target]
        merges[tip] = target
    return merges

def getNetworkX(mycp, geometry = None, hyphae_only = True):
    # undirected networkx graph of the hyphal network (inside geometry, e.g. small_hyphae_dish)
    # nodes: node ids of the plant (cut nodes from crop get new ids), attribute "pos" (x, y, z)
    # edges: one per segment, attributes "length" (3D), "length_xy" (projected to x-y plane) and all per segment data of getMycSegmentAnalyser
    # anastomoses: the tip node is merged with the node it fused with, so the hyphal trees are connected
    import networkx as nx
    ana = getMycSegmentAnalyser(mycp)
    if hyphae_only:
        ana.filter("organType", 5)
    if geometry is not None:
        ana.crop(geometry) # no pack(), it would renumber the nodes and the anastomosis ids would not fit anymore
    nodes = np.array([[n.x, n.y, n.z] for n in ana.nodes])
    segs = np.array([[s.x, s.y] for s in ana.segments], dtype=int)
    inside = set(segs.ravel().tolist())
    merges = {tip: target for tip, target in getAnastomosisMerges(mycp).items() if target in inside} # fusions with hyphae outside geometry are not connected
    segs_merged = np.vectorize(lambda i: merges.get(i, i), otypes=[int])(segs) if (merges and len(segs) > 0) else segs
    data = {"length": np.array(ana.getParameter("length"))}
    data["length_xy"] = np.linalg.norm(nodes[segs[:, 0], :2] - nodes[segs[:, 1], :2], axis=1) if len(segs) > 0 else np.zeros(0)
    for key in ana.data.keys():
        data[key] = np.array(ana.data[key])
    G = nx.Graph()
    for i in np.unique(segs_merged):
        G.add_node(int(i), pos = tuple(nodes[i].tolist()))
    keys = list(data.keys())
    for k, (u, v) in enumerate(segs_merged):
        if u != v: # an anastomosis segment of zero length could become a self loop
            G.add_edge(int(u), int(v), **{key: float(data[key][k]) for key in keys})
    return G

def lengthWeight(dim):
    # edge attribute used as weight: 3D length or length projected to the x-y plane (dish seen from above)
    return "length" if dim == 3 else "length_xy"

def simplifyGraph(G, weights = ("length", "length_xy")):
    # multigraph with only junctions and ends of G (degree != 2), each chain of degree 2 nodes becomes one edge
    # with summed weights and the number of segments ("segments"); a multigraph keeps parallel chains and loops,
    # so total length is unchanged and shortest paths between the kept nodes are the same as in G
    import networkx as nx
    keep = {n for n, d in G.degree() if d != 2}
    for comp in nx.connected_components(G): # closed loops without any junction
        if not any(n in keep for n in comp):
            keep.add(min(comp))
    S = nx.MultiGraph()
    for n in keep:
        S.add_node(n, **G.nodes[n])
    visited = set()
    for n in keep:
        for nb in G.neighbors(n):
            if frozenset((n, nb)) in visited:
                continue
            sums = {w: 0. for w in weights}
            segments = 0
            prev, cur = n, nb
            while True:
                visited.add(frozenset((prev, cur)))
                for w in weights:
                    sums[w] += G.edges[prev, cur][w]
                segments += 1
                if cur in keep:
                    break
                prev, cur = cur, next(x for x in G.neighbors(cur) if x != prev) # cur has degree 2
            S.add_edge(n, cur, segments = segments, **sums)
    return S

def getReferenceNetworks(S, dim = 3):
    # Delaunay triangulation and Euclidean minimum spanning tree on the nodes of S (use simplifyGraph first),
    # dim = 2: triangulation of the x-y positions (dish seen from above), dim = 3: tetrahedra in x, y, z
    # edges get "length" and "length_xy" (straight lines), the MST minimises lengthWeight(dim)
    import networkx as nx
    from scipy.spatial import Delaunay
    nodes = list(S.nodes)
    pos = np.array([S.nodes[n]["pos"] for n in nodes])
    tri = Delaunay(pos[:, :dim], qhull_options = "QJ") # QJ: robust against duplicate and coplanar points
    D = nx.Graph()
    for n in nodes:
        D.add_node(n, **S.nodes[n])
    for simplex in tri.simplices:
        for a in range(len(simplex)):
            for b in range(a + 1, len(simplex)):
                i, j = simplex[a], simplex[b]
                D.add_edge(nodes[i], nodes[j], length = float(np.linalg.norm(pos[i] - pos[j])), length_xy = float(np.linalg.norm(pos[i, :2] - pos[j, :2])))
    M = nx.minimum_spanning_tree(D, weight = lengthWeight(dim)) # the Euclidean MST is always part of the Delaunay graph
    return D, M

def weightedEfficiency(G, dim = 3):
    # global efficiency of Latora & Marchiori (2001): sum of 1/(shortest path length) over all node pairs,
    # normalised by sum of 1/(straight line distance), 1 = every pair connected by a straight line;
    # unconnected pairs contribute 0, pairs at the same position are skipped
    import networkx as nx
    from scipy.spatial.distance import pdist
    weight = lengthWeight(dim)
    nodes = list(G.nodes)
    pos = np.array([G.nodes[n]["pos"] for n in nodes])[:, :dim]
    eucl = pdist(pos)
    inv_eucl = 2 * np.sum(1. / eucl[eucl > 0])
    inv_path = 0.
    for s in nodes:
        dist = nx.single_source_dijkstra_path_length(G, s, weight = weight)
        inv_path += sum(1. / d for t, d in dist.items() if t != s and d > 0)
    return float(inv_path / inv_eucl)

def getSourcesTargets(G, barrier_thickness, tol = 1.e-5):
    # sources: nodes where hyphae cross into the hyphal compartment, i.e. the nodes crop(small_hyphae_dish) created
    # on the cut plane x = -barrier_thickness/2 (degree 1); targets: all other nodes of degree 1 (tips and ends)
    ends = [n for n, d in G.degree() if d == 1]
    sources = [n for n in ends if abs(G.nodes[n]["pos"][0] + barrier_thickness / 2) < tol]
    s = set(sources)
    targets = [n for n in ends if n not in s]
    return sources, targets

def weightedEfficiencySubset(G, sources, targets = None, dim = 3):
    # efficiency of Latora & Marchiori (2001) only for the pairs (source, target), or for all pairs within sources if targets is None:
    # sum of 1/(shortest path length) normalised by sum of 1/(straight line distance), unconnected pairs contribute 0
    import networkx as nx
    from scipy.spatial.distance import cdist
    weight = lengthWeight(dim)
    sources = list(sources)
    if targets is None:
        targets = sources
    else:
        targets = [t for t in targets if t not in set(sources)]
    pos_s = np.array([G.nodes[n]["pos"] for n in sources])[:, :dim]
    pos_t = np.array([G.nodes[n]["pos"] for n in targets])[:, :dim]
    eucl = cdist(pos_s, pos_t) # pairs within one set are counted twice, as are the paths below
    inv_eucl = np.sum(1. / eucl[eucl > 0])
    inv_path = 0.
    for s in sources:
        dist = nx.single_source_dijkstra_path_length(G, s, weight = weight)
        inv_path += sum(1. / dist[t] for t in targets if t != s and t in dist and dist[t] > 0)
    return float(inv_path / inv_eucl)

def getTerminalRootNodes(mycp, S, barrier_thickness, root_dist = 0.5):
    # node sets as in Oyarte Galvez et al. (2025), both restricted to the nodes of S (e.g. simplifyGraph(getNetworkX(mycp, small_hyphae_dish))):
    # R (root nodes): all nodes within root_dist [cm] (5 mm) of the barrier face on the hyphal side (x = barrier_thickness/2),
    #                 i.e. sites of crossing from the root compartment
    # T (terminal nodes): growing tips (last node of active non-BAS hyphae) and attachment sites of BAS (first node of a BAS,
    #                 subType 3, whose parent is no BAS); nodes in R are not in T
    x_root = barrier_thickness / 2 + root_dist
    R = [n for n in S.nodes if S.nodes[n]["pos"][0] <= x_root]
    in_R = set(R)
    T = set()
    for h in mycp.getOrgans(pb.hyphae):
        if int(h.getParameter("subType")) == 3:
            p = h.getParent()
            if not (int(p.getParameter("organType")) == 5 and int(p.getParameter("subType")) == 3):
                T.add(h.getNodeId(0))
        elif h.isActive():
            T.add(h.getNodeId(h.getNumberOfNodes() - 1))
    T = [n for n in T if n in S and n not in in_R]
    return T, R

def getEfficiencies(S, T, R, dim = 3):
    # (i) global efficiency: transport between arbitrary pairs of the transport relevant nodes T and R
    # (ii) root efficiency: transport between the terminal nodes T and the root nodes R
    return weightedEfficiencySubset(S, list(T) + list(R), None, dim), weightedEfficiencySubset(S, T, R, dim)

def betweennessSubset(G, sources, targets, dim = 3, edges = False):
    # fraction of all shortest (source, target) paths running through each node (or edge, if edges = True);
    # networkx counts every path in an undirected graph with 0.5 and its normalized=True divides by (n-1)(n-2), so rescale here
    import networkx as nx
    weight = lengthWeight(dim)
    if edges:
        bc = nx.edge_betweenness_centrality_subset(G, sources, targets, normalized = False, weight = weight)
    else:
        bc = nx.betweenness_centrality_subset(G, sources, targets, normalized = False, weight = weight)
    scale = 2. / (len(sources) * len(targets)) if (sources and targets) else 0.
    return {k: v * scale for k, v in bc.items()}