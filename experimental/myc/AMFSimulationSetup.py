import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import plantbox as pb
import time

import AMFAnalysis as amf
import medium_heights as mh

# Geometry of the split petri dish and the simulation of roots and hyphae in it (root phase, waiting for t*, hyphae phase)

def defaultPetridish(nRings = 15, hyphae_layer = 0.1, opening_at_top = True, rootradius = 0.1, verbose = True):
    # split petri dish of the experiment (Oyarte Galvez et al. 2025): 94 mm dish, barrier 1.6 mm, opening 50 x 2 mm;
    # the modelled domain is the medium (z in [-height, 0]), its height follows from the filled volumes of the experiment (medium_heights.py)
    # hyphae_layer: hyphae grow on the cellophane in the fungal compartment, i.e. in a thin layer at the top (None = whole height)
    # returns a dict with the parameters and the geometries of makedishes ("petri_dish", "small_hyphae_dish", "half_dish", "rings")
    diameter = 9.4
    barrier_thickness = 0.16
    opening_length = 5.0
    opening_height = 0.2
    heights = mh.estimate_heights(diameter, barrier_thickness)
    if verbose:
        mh.print_heights(heights, diameter, barrier_thickness)
    height = heights["root"]
    barrier_height = height
    petri_dish, small_hyphae_dish, half_dish, rings = makedishes(diameter, height, barrier_thickness, barrier_height, opening_length, opening_height,
                                                                 rootradius, nRings, opening_at_top = opening_at_top, hyphae_layer = hyphae_layer)
    return {"diameter": diameter, "radius": diameter / 2, "height": height, "barrier_thickness": barrier_thickness, "barrier_height": barrier_height,
            "opening_length": opening_length, "opening_height": opening_height, "hyphae_layer": hyphae_layer, "nRings": nRings,
            "petri_dish": petri_dish, "small_hyphae_dish": small_hyphae_dish, "half_dish": half_dish, "rings": rings}

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
    # (a hypha whose ancestor died stays active but is not simulated anymore, see AMFDensities.getActiveTipsPerRing)
    t_min = plant.getSimTime() - window
    n = 0
    for o in plant.getOrgans(ot):
        k = o.getNumberOfNodes() - 1
        if o.getNodeCT(k) >= t_min and geometry.getDist(o.getNode(k)) < 0:
            n += 1
    return n

def makeSimulation(seed, path, name, height, petri_dish, small_hyphae_dish, half_dish, dt, nSteps, hours_hyphae, filename, animation = False, verbose = False,
                   max_wait_days = None, stats = None, min_growing_tips = None, stall_fraction = None, growing_window = 0.25, animation_dir = "animation/",
                   turn_off_root_growth = False, turn_off_colonization = False, progress_every = None,
                   seed_pos_x = -9.4 / 2 + 0.5, hyphal_emergence_density = 4.0, lmbd = 0.15, crossing_length = 10.):
    # set up simulation
    # max_wait_days: stop waiting for t* (crossing_length cm hyphae behind the barrier) after this many simulated days, then returns mycp, None, None
    #                (None = wait without limit, as before). A limit in simulated time is reproducible, a wall clock limit is not.
    # min_growing_tips: stop the hyphae phase after t* as soon as fewer tips than this grow behind the barrier (created a node in small_hyphae_dish
    #                   within the last growing_window days), then returns mycp, None, None with status "stalled" (None = no check)
    # stall_fraction: same, if the number of growing tips falls below stall_fraction times its maximum so far (e.g. 0.5, None = no check)
    # turn_off_root_growth, turn_off_colonization: at t*, when the hyphae have grown to the other side of the barrier, the roots stop growing,
    #                   or stop colonization and the emergence of new hyphae from the roots (mycp.turnOffRoots), the hyphae keep growing
    # progress_every: print one line with the state of the simulation every progress_every time steps (None = no output), e.g. 24 = once per simulated day
    # seed_pos_x: x position of the seed [cm] (default 5 mm from the edge of the 94 mm dish), it lies at y = 0, z = -height/2
    # hyphal_emergence_density, lmbd: set for all root types after the root phase, when the hyphae start to emerge
    # crossing_length: hyphal length [cm] behind the barrier that defines t*
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
    seed_para.seedPos.x = seed_pos_x
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
            amf.getMycSegmentAnalyser(mycp,filename = animation_dir + filename + "_anim",step = i, stdwrite= True)    # look at roots and container
        progress("root phase", i)

    root = mycp.getOrganRandomParameter(pb.root)
    for rp in root:
        rp.hyphalEmergenceDensity = hyphal_emergence_density
        rp.lmbd = lmbd
        mycp.setOrganRandomParameter(rp)

    mycp.changeGeometry(5,petri_dish)

    wall_root = time.perf_counter() - t_wall
    t_wall = time.perf_counter()
    t_root_end = mycp.getSimTime()
    crossed_t, crossed_len = [], []
    status = "ok"
    crossed_barrier = 0

    step = 0
    while crossed_barrier < crossing_length:
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
        progress("waiting for t*", step, f"behind barrier {crossed_barrier:.2f} of {crossing_length:g} cm")
        step += 1

    if stats is not None:
        stats.update({"status": status, "wall_root": wall_root, "wall_wait": time.perf_counter() - t_wall, "wall_hyphae": 0.,
                      "t_root_end": t_root_end, "wait_days": mycp.getSimTime() - t_root_end,
                      "crossed_t": np.array(crossed_t), "crossed_len": np.array(crossed_len)})
    if status != "ok":
        print("Seed " + str(seed) + ": less than " + f"{crossing_length:g}" + " cm hyphae behind the barrier after " + str(max_wait_days) + " days, stopped")
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
            amf.getMycSegmentAnalyser(mycp,filename = animation_dir + filename + "_anim",step = i+nSteps, stdwrite= True)    # look at roots and container
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
