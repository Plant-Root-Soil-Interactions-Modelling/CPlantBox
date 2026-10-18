import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import plantbox as pb
import time
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import plantbox.visualisation.vtk_plot as vp

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

def getParaSumperRing(parameter, times, plant, rings):
    print("Parameter sum for "+ parameter)
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

def makeBoxes(radius,nRings,height, other = 0.):
    length = radius*np.sqrt(1/nRings)
    box = pb.SDF_PlantContainer(length,length,height,True)

def setUpSimulationTime(simTime, fps):
    # set up simulation time
    simTime = simTime
    fps = fps
    dt = 1 / fps
    nSteps = int(simTime / dt)
    return dt, nSteps

def makeSimulation(seed, path, name, height, petri_dish, small_hyphae_dish, half_dish, dt, nSteps, hours_hyphae, filename, animation = False, verbose = False):
    # set up simulation
    #### WORK IN PROGRESS
    mycp = pb.MycorrhizalPlant(seed)
    mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)
    seed_para = pb.SeedRandomParameter(mycp)
    seed_para.seedPos.z = -height / 2
    seed_para.seedPos.x = -9.4 / 2 + 0.5
    seed_para.seedPos.y = 0
    mycp.setOrganRandomParameter(seed_para)
    mycp.setGeometry(half_dish)
    mycp.initialize()

    print("Starting simulation with seed: " + str(seed))
    for i in range(nSteps):
        if (i % 10 == 0):
            print("Step " + str(i) + " of " + str(nSteps))
        mycp.simulate(dt,verbose)
        if (animation):
            getMycSegmentAnalyser(mycp,filename = "animation/" + filename + "_anim",step = i, stdwrite= True)    # look at roots and container

    root = mycp.getOrganRandomParameter(pb.root)
    for rp in root:
        rp.hyphalEmergenceDensity = 4.0
        rp.lmbd = 0.15
        mycp.setOrganRandomParameter(rp)

    mycp.changeGeometry(5,petri_dish)

    crossed_barrier = 0

    while crossed_barrier < 10:
        # N+=1
        mycp.simulate(dt,verbose)
        ana = getMycSegmentAnalyser(mycp)
        ana.filter("organType",5)
        ana.pack()
        crossed_barrier = ana.getSummed("length",small_hyphae_dish)

    time_star = mycp.getSimTime()
    for i in range(0, hours_hyphae):
        if (i % 10 == 0):
            print("Step " + str(i) + " of " + str(hours_hyphae))
        mycp.simulate(dt,verbose)
        if (animation):
            getMycSegmentAnalyser(mycp,filename = "animation/" + filename + "_anim",step = i+nSteps, stdwrite= True)    # look at roots and container
    times = np.linspace(time_star, max(mycp.getParameter("creationTime"))+0.01, 100)
    # times = times - time_star
    return mycp, time_star, times

def makeFreeSimulation(seed, path, name, dt, nSteps, filename, animation = False, verbose = False):
    # set up simulation
    #### WORK IN PROGRESS
    mycp = pb.MycorrhizalPlant(seed)
    mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)
    seed_para = pb.SeedRandomParameter(mycp)
    seed_para.seedPos.z = - 3
    mycp.setOrganRandomParameter(seed_para)
    root = mycp.getOrganRandomParameter(pb.root)
    for rp in root:
        rp.hyphalEmergenceDensity = 4.0
        rp.lmbd = 0.15
        mycp.setOrganRandomParameter(rp)
    mycp.initialize()

    print("Starting simulation with seed: " + str(seed))
    for i in range(nSteps):
        if (i % 10 == 0):
            print("Step " + str(i) + " of " + str(nSteps))
        mycp.simulate(dt,verbose)
        if (animation):
            getMycSegmentAnalyser(mycp,filename = "animation/" + filename + "_anim" +str(i) + ".vtp", stdwrite= True)    # look at roots and container
    times = np.linspace(0, max(mycp.getParameter("creationTime"))+0.01, 100)
    return mycp, times

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

def findAnaDensDep(t, ana, rings, power = 1):
    ana_densities = getParaSumperRing("anastomosis", t, ana, rings)
    tip_densities = getParaSumperRing("nodeTips", t, ana, rings)
    length_densities = getParaSumperRing("length", t, ana, rings)
    ana_rates = np.zeros(len(rings),len(t[1:]))
    for i in range(len(rings)):
        ana_rates[i,:] = ana_densities[i,:]/tip_densities[i,:]
        coefs = np.polyfit(length_densities[i,:], ana_rates[i,:], deg = power)

def prepare_data(lenMat, vol, anaMat, tipMat):
    rho = lenMat / vol
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

# def sensitivity_analysis(seed,path,name,runnerhyphae = {},BAS = {}, geometry = {},run = {}):
#     mycp = pb.MycorrhizalPlant(seed)
#     mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)
#     seed_para = pb.SeedRandomParameter(mycp)
#     seed_para.seedPos.z = -geometry["height"] / 2
#     seed_para.seedPos.x = -9.4 / 2 + 0.5
#     seed_para.seedPos.y = 0
#     mycp.setOrganRandomParameter(seed_para)
#     hyphae = mycp.getOrganRandomParameter(pb.hyphae)
#     for h in hyphae:
#         if h.subType == 1:
#             h.a = runnerhyphae["a"]
#         if h.subType == 2:
#             h.a = BAS["a"]

#     mycp.setGeometry(geometry["half_dish"])
#     mycp.initialize()

#     print("Starting simulation with seed: " + str(seed))
#     for i in range(run["nSteps"]):
#         if (i % 10 == 0):
#             print("Step " + str(i) + " of " + str(run["nSteps"]))
#         mycp.simulate(run["dt"],run["verbose"])
#     root = mycp.getOrganRandomParameter(pb.root)
#     for rp in root:
#         rp.hyphalEmergenceDensity = 4.0
#         rp.lmbd = 0.15
#         mycp.setOrganRandomParameter(rp)

#     mycp.changeGeometry(5,geometry["petri_dish"])

#     crossed_barrier = 0

#     while crossed_barrier < 3:
#         # N+=1
#         mycp.simulate(dt,verbose)
#         ana = getMycSegmentAnalyser(mycp)
#         ana.filter("organType",5)
#         ana.filter("subType",1,2)
#         ana.pack()
#         crossed_barrier = ana.getSummed("nodeTips",small_hyphae_dish)

#     crossed_time = mycp.getSimTime()
#     for i in range(0, hours_hyphae):
#         if (i % 10 == 0):
#             print("Step " + str(i) + " of " + str(hours_hyphae))
#         mycp.simulate(dt,verbose)
#         if (animation):
#             getMycSegmentAnalyser(mycp,filename = "animation/" + filename + "_anim" +str(i) + ".vtp", stdwrite= True)    # look at roots and container
#     times = np.linspace(crossed_time, max(mycp.getParameter("creationTime"))+0.01, 100)
#     return mycp, crossed_time, times

# def estimate_carbon_cost(length,radius,pdry,pcarbon,dcell,CUE,)