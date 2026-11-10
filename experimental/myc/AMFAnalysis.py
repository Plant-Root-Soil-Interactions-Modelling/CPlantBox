import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import plantbox as pb
import time
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import networkx as nx
import glob
import os
import warnings
from scipy.interpolate import interp1d
from scipy.spatial import Delaunay
from scipy.spatial.distance import pdist
from scipy.spatial.distance import cdist

def petridishOutputDirs():
    # output folders of the petri dish scripts (petridishsimauto.py, petridishcombine.py, petridishtiming.py, petridishstats.py):
    # below the environment variable PETRIDISH_OUT (set by run_petridish.sh on the cluster), else below petridish_results in the
    # current folder. Returns a dict with "base" (combined results), "densities" (seed_*.npz), "plots" and "animation" (vtp), each ending with "/"
    base = os.environ.get("PETRIDISH_OUT", "petridish_results")
    dirs = {"base": base, "densities": os.path.join(base, "densities"), "plots": os.path.join(base, "plots"), "animation": os.path.join(base, "animation")}
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)
    return {k: os.path.join(d, "") for k, d in dirs.items()}

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

HYPHAL_SUBTYPES = {1: "runner", 2: "bas"} # hyphal sub types of TwoHyphaePlusBAS.xml (RunnerHyphae, BAS)

def getParasPerRing(parameters, times, plant, rings, subTypes = None, increments = ()):
    # same matrices as getParaSumperRing (rings x len(times)-1), but for all parameters at once:
    # one SegmentAnalyser, ring of each segment computed once, time intervals binned with numpy
    # subTypes: hyphal sub types (e.g. [1, 2]), additionally returns parameter + "_st" + subType with only the hyphae of that sub type
    # increments: parameters for which additionally parameter + "_new" (and parameter + "_st" + subType + "_new") is returned:
    #             the amount created in each interval (times[i], times[i+1]] instead of the cumulative amount
    ana = getMycSegmentAnalyser(plant)
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
                mat[:, 0] = summed[:, 1] # first column only counts (times[0], times[1]], as in getParaSumperRing
            mats[parameter + suffix] = mat
            if parameter in increments:
                mats[parameter + suffix + "_new"] = summed[:, 1:]
    return mats

def getParaSumperRing(parameter, times, plant, rings):
    return getParasPerRing([parameter], times, plant, rings)[parameter]

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

def makedishes(diameter, height, barrier_thickness, barrier_height, opening_length, opening_height, rootradius,nRings,
               opening_at_top = False, hyphae_layer = None):
    # height: height of the modelled domain, z in [-height, 0] (for the medium heights of the experiment see medium_heights.py)
    # opening_at_top: opening at the top of the barrier, z in [-opening_height, 0], instead of at the bottom
    # hyphae_layer: thickness of the layer at the top in which the hyphae grow in the fungal compartment (on the cellophane),
    #               None = whole height; small_hyphae_dish and rings are restricted to this layer
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
    if opening_at_top:
        bottompoint = -bottompoint # opening at the upper edge of the barrier
    # container for opening moved to the right position
    helper_staff2 = pb.SDF_RotateTranslate(helper_staff2, 0,0,pb.Vector3d(0, bottompoint, midpoint))
    # opening made in the barrier
    helper_dish2 = pb.SDF_Difference(helper_staff, helper_staff2)
    # barrier moved to the right position
    moved_helper_dish2 = pb.SDF_RotateTranslate(helper_dish2, 90, pb.SDF_Axis.xaxis , pb.Vector3d(0, -radius, -height/2))
    petri_dish = pb.SDF_Difference(petri_dish, moved_helper_dish2)
    h_hyphae = height
    if hyphae_layer is not None:
        # behind the barrier (x > barrier_thickness/2) the hyphae can only grow in the top layer, z in [-hyphae_layer, 0]
        h_hyphae = hyphae_layer
        below_layer = pb.SDF_PlantContainer(radius, radius, height - hyphae_layer, True)
        below_layer = pb.SDF_RotateTranslate(below_layer, 0, 0, pb.Vector3d(radius + barrier_thickness/2, 0, -hyphae_layer))
        petri_dish = pb.SDF_Difference(petri_dish, below_layer)

    moved_helper_dish_hyphae = pb.SDF_RotateTranslate(helper_dish, 0, 0, pb.Vector3d(-(radius+barrier_thickness/2), 0, 0))
    small_dish = pb.SDF_PlantContainer(radius,radius,h_hyphae,False)
    small_hyphae_dish = pb.SDF_Difference(small_dish, moved_helper_dish_hyphae)

    small_dish = pb.SDF_PlantContainer(radius*np.sqrt(1/nRings),radius*np.sqrt(1/nRings),h_hyphae,False)
    ringone = pb.SDF_Difference(small_dish, moved_helper_dish_hyphae)
    rings = []
    rings.append(ringone)
    for i in range(2, nRings+1):
        small_dish = pb.SDF_PlantContainer(radius*np.sqrt(i/nRings),radius*np.sqrt(i/nRings),h_hyphae,False)
        small_hyphae_dish2 = pb.SDF_Difference(small_dish, moved_helper_dish_hyphae)
        old_dish = pb.SDF_Difference(pb.SDF_PlantContainer(radius*np.sqrt((i-1)/nRings),radius*np.sqrt((i-1) /nRings),h_hyphae,False),moved_helper_dish_hyphae)
        small_hyphae_dish2 = pb.SDF_Difference(small_hyphae_dish2,old_dish)
        rings.append(small_hyphae_dish2)
    return petri_dish, small_hyphae_dish, half_dish, rings

def makeBoxes(radius, height, boxlength, barrier_thickness, small_hyphae_dish, resolution = 0.005):
    # non overlapping square boxes (edge length boxlength) tiling the hyphae side of the dish (x > -barrier_thickness/2),
    # each box is cut with small_hyphae_dish, boxes outside the dish are dropped
    # returns boxes, centres of the uncut boxes (n x 2) and volume of each cut box
    # height: height of small_hyphae_dish, i.e. hyphae_layer if makedishes was called with it
    box =pb.SDF_PlantContainer(boxlength/2, boxlength/2, height, True) # square container has half edge length r
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
    # height: height of the rings, i.e. hyphae_layer if makedishes was called with it; height = 1 gives the area of each ring
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

def countGrowingTips(plant, geometry, window, ot = 5):
    # number of organs whose last node lies in geometry and was created within the last window days, i.e. tips that are really growing
    # (a hypha whose ancestor died stays active but is not simulated anymore, see getActiveTipsPerRing)
    t_min = plant.getSimTime() - window
    n = 0
    for o in plant.getOrgans(ot):
        k = o.getNumberOfNodes() - 1
        if o.getNodeCT(k) >= t_min and geometry.getDist(o.getNode(k)) < 0:
            n += 1
    return n

def makeSimulation(seed, path, name, height, petri_dish, small_hyphae_dish, half_dish, dt, nSteps, hours_hyphae, filename, animation = False, verbose = False,
                   max_wait_days = None, stats = None, min_growing_tips = None, stall_fraction = None, growing_window = 0.25, animation_dir = "animation/",
                   turn_off_root_growth = False, turn_off_colonization = False, progress_every = None):
    # set up simulation
    # max_wait_days: stop waiting for t* (10 cm hyphae behind the barrier) after this many simulated days, then returns mycp, None, None
    #                (None = wait without limit, as before). A limit in simulated time is reproducible, a wall clock limit is not.
    # min_growing_tips: stop the hyphae phase after t* as soon as fewer tips than this grow behind the barrier (created a node in small_hyphae_dish
    #                   within the last growing_window days), then returns mycp, None, None with status "stalled" (None = no check)
    # stall_fraction: same, if the number of growing tips falls below stall_fraction times its maximum so far (e.g. 0.5, None = no check)
    # turn_off_root_growth, turn_off_colonization: at t*, when the hyphae have grown to the other side of the barrier, the roots stop growing,
    #                   or stop colonization and the emergence of new hyphae from the roots (mycp.turnOffRoots), the hyphae keep growing
    # progress_every: print one line with the state of the simulation every progress_every time steps (None = no output), e.g. 24 = once per simulated day
    # stats: optional dict, filled with "status" ("ok", "not_crossed" or "stalled"), the wall clock time of each phase [s] ("wall_root", "wall_wait",
    #        "wall_hyphae"), "t_root_end" (simulation time after the root phase), "wait_days" (simulated days until t*, or until the stop),
    #        the length behind the barrier after each waiting step ("crossed_t", "crossed_len") and, with min_growing_tips, the number of
    #        growing tips behind the barrier after each hyphae step ("growing_t", "growing_n")
    #### WORK IN PROGRESS
    t_wall = time.perf_counter()
    t_start = t_wall
    def progress(phase, i, extra = ""):
        # one line with the state after time step i of the phase, flushed so it shows up immediately in slurm output files
        if progress_every is None or (i + 1) % progress_every != 0:
            return
        print(f"seed {seed} | {phase} step {i + 1} | t {mycp.getSimTime():.2f} d | wall {time.perf_counter() - t_start:.0f} s | "
              f"roots {mycp.getSummed('length', pb.root):.1f} cm | hyphae {len(mycp.getOrgans(pb.hyphae))}, {mycp.getSummed('length', pb.hyphae):.1f} cm"
              + (" | " + extra if extra else ""), flush=True)
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
            getMycSegmentAnalyser(mycp,filename = animation_dir + filename + "_anim",step = i, stdwrite= True)    # look at roots and container
        progress("root phase", i)

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

    step = 0
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
        progress("waiting for t*", step, f"behind barrier {crossed_barrier:.2f} of 10 cm")
        step += 1

    if stats is not None:
        stats.update({"status": status, "wall_root": wall_root, "wall_wait": time.perf_counter() - t_wall, "wall_hyphae": 0.,
                      "t_root_end": t_root_end, "wait_days": mycp.getSimTime() - t_root_end,
                      "crossed_t": np.array(crossed_t), "crossed_len": np.array(crossed_len)})
    if status != "ok":
        print("Seed " + str(seed) + ": less than 10 cm hyphae behind the barrier after " + str(max_wait_days) + " days, stopped")
        return mycp, None, None
    t_wall = time.perf_counter()

    time_star = mycp.getSimTime()
    if progress_every is not None:
        print(f"seed {seed} | t* reached at t {time_star:.2f} d, {time_star - t_root_end:.2f} d after the root phase", flush=True)
    if turn_off_root_growth or turn_off_colonization:
        mycp.turnOffRoots(turn_off_root_growth, turn_off_colonization)
    growing_t, growing_n = [], []
    for i in range(0, hours_hyphae):
        # if (i % 10 == 0):
        #     print("Step " + str(i) + " of " + str(hours_hyphae))
        mycp.simulate(dt,verbose)
        if (animation):
            getMycSegmentAnalyser(mycp,filename = animation_dir + filename + "_anim",step = i+nSteps, stdwrite= True)    # look at roots and container
        if min_growing_tips is not None or stall_fraction is not None:
            growing_t.append(mycp.getSimTime())
            growing_n.append(countGrowingTips(mycp, small_hyphae_dish, growing_window))
            if (min_growing_tips is not None and growing_n[-1] < min_growing_tips) or \
               (stall_fraction is not None and growing_n[-1] < stall_fraction * max(growing_n)):
                status = "stalled"
                break
        progress("hyphae after t*", i, f"{mycp.getSimTime() - time_star:.2f} d after t*"
                 + (f" | growing tips behind barrier {growing_n[-1]} (max {max(growing_n)})" if growing_n else ""))
    if stats is not None:
        stats.update({"status": status, "wall_hyphae": time.perf_counter() - t_wall, "growing_t": np.array(growing_t), "growing_n": np.array(growing_n)})
    if status != "ok":
        print("Seed " + str(seed) + ": only " + str(growing_n[-1]) + " growing tips behind the barrier " + str(mycp.getSimTime() - time_star) + " days after t*, stopped")
        return mycp, None, None
    nmbtimes = (max(mycp.getParameter("creationTime")) - time_star)*24
    nmbtimespoint = (max(mycp.getParameter("creationTime")) + 0.1 - time_star)*24
    times = np.linspace(time_star, max(mycp.getParameter("creationTime"))+0.01, max(round(nmbtimes),round(nmbtimespoint)))
    # times = times - time_star
    return mycp, time_star, times



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


def interpolate_to_common_time(simulations, key, n_times=99):
    # interpolates sim[key] (rings x times) of all simulations onto the time range they share

    t_start = max(sim["times"][0] for sim in simulations)
    t_end = min(sim["times"][-1] for sim in simulations)
    if t_start >= t_end:
        raise ValueError(f"interpolate_to_common_time: time ranges of the simulations do not overlap (latest start {t_start}, earliest end {t_end})")
    t_common = np.linspace(t_start, t_end, n_times)
    interp = np.stack([interp1d(sim["times"], sim[key], axis=1)(t_common) for sim in simulations])
    return t_common, interp   # (n_sim, n_r, n_t)


def save_simulation(filename, sim):
    # stores one simulation dict (only arrays/numbers, e.g. "times", "tip_dens", ...) as .npz
    np.savez(filename, **sim)

def load_simulations(pattern):
    # loads all .npz files matching pattern (e.g. "results/seed_*.npz") into a list of dicts, sorted by filename

    simulations = []
    for f in sorted(glob.glob(pattern)):
        with np.load(f) as data:
            simulations.append({key: data[key] for key in data.files})
    # print("Loaded " + str(len(simulations)) + " simulations from " + pattern)
    return simulations


def load_petridish_results(pattern, verbose=True):
    # loads the per-seed results of petridishsimauto.py (e.g. "petridish_results/densities/seed_*.npz") that can be combined:
    # seeds stopped before t* (status "not_crossed", see max_wait_days) or because the front stalled after t* (status "stalled",
    # see min_growing_tips, stall_fraction) have no densities and are skipped, older files have no status.
    # Older files have densities per cm^3: with hyphae_layer (thin layer geometry) they are converted to mm^-2
    # (cm^-3 * layer [cm] = cm^-2 = mm^-2 * 100), without it (old geometry, densities per volume of the whole height) they are skipped.
    # "times" are shifted to days since t* (the absolute times of different seeds do not overlap)
    simulations = load_simulations(pattern)
    stopped = [int(sim["seed"]) for sim in simulations if str(sim.get("status", "ok")) != "ok"]
    simulations = [sim for sim in simulations if str(sim.get("status", "ok")) == "ok"]
    for sim in simulations:
        if "density_per" not in sim and "hyphae_layer" in sim:
            h = float(sim["hyphae_layer"])
            for key in ["tip_dens", "stopped_tip_dens"]:
                sim[key] = sim[key] * h / 100
            sim["length_densities"] = sim["length_densities"] * h / 10 # cm cm^-3 -> cm cm^-2 -> mm mm^-2
            sim["density_per"] = "mm^2"
    old = [int(sim["seed"]) for sim in simulations if "density_per" not in sim]
    simulations = [sim for sim in simulations if "density_per" in sim]
    for sim in simulations:
        sim["times"] = sim["times"] - sim["time_star"]
    if verbose:
        if old:
            print("skipping " + str(len(old)) + " seeds of the old geometry (rerun them): " + str(old))
        print(str(len(simulations)) + " simulations reached t*, " + str(len(stopped)) + " were stopped: " + str(stopped))
    return simulations

def estimate_carbon_cost_per_Org(length,radius,pdry,pcarbon,dcell,CUE):
    gamma = np.pi * radius *radius*dcell*pdry*pcarbon
    gamma = gamma / CUE
    carbon_cost = gamma*length
    return carbon_cost
