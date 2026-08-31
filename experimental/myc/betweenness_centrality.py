import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.visualisation import figure_style
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
import time
import matplotlib as mpl
import AMFAnalysis as amf

# I changed the c++ code, so set might not be valid anymore
good_seeds = [10, 20, 30, 60, 70, 80, 90, 100] # these are seeds where hyphae cross the barrier, but have not been assessed in any other way for any required behaviour


def makesimulation(seed):
    mycp = pb.MycorrhizalPlant(seed)
    path = "tomatoparameters/"
    name = "TwoHyphaePlusBAS"

    start = time.perf_counter()
    animation = False
    mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)

    ### initial root parameters
    root = mycp.getOrganRandomParameter(pb.root)
    for rp in root:
        rp.dx = 0.1
        rp.maxAge = 100 # maximal colonization age
        # rp.a = 0.01
        mycp.setOrganRandomParameter(rp)

    ## Setting up petri dish
    diameter = 9.4
    height = 1.6
    # introduce parameters for barrier and opening
    barrier_thickness = 0.16
    barrier_height = height
    opening_length = 5.0
    opening_height = 0.2

    nRings = 25
    petri_dish, small_hyphae_dish, half_dish, rings = amf.makedishes(diameter, height, barrier_thickness, barrier_height, opening_length, opening_height, root[0].a,nRings)


    # make sure to set the seed position to 0 because of the petri dish
    seed_parameter = pb.SeedRandomParameter(mycp)
    seed_parameter.seedPos.z = -height / 6 # seed is positioned in the middle of the petri dish in the z direction
    seed_parameter.seedPos.x = - 1.0
    seed_parameter.seedPos.y = 0
    mycp.setOrganRandomParameter(seed_parameter)

    mycp.setGeometry(half_dish)
    mycp.initialize(True)

    # set up simulation times etc.
    simtime = 5
    fps = 24
    N = fps * simtime
    dt = simtime / N

    # filename = "splitpetri_dish_parametrisation_" + str(simtime)

    filename = "petridish_seed_" + str(seed)

    print("Starting simulation with seed " + str(seed) + " and filename " + filename)
    # Start simulation
    start = time.perf_counter()
    for i in range(0, N):
        if i % 6 == 0:
            print("Step " + str(i) + " of " + str(N))
        mycp.simulate(dt,True)
        if (animation):
            ana = amf.getMycSegmentAnalyser(mycp)
            ana.write("animation/" + filename + "_hoursBCB_" +str(i) + ".vtp", ["radius", "subType", "creationTime", "organType", "colonization", "colonizationTime", "anastomosis"])
    # look at roots and container
    # vp.plot_roots_and_container(mycp,half_dish)


    afterroots = time.perf_counter()
    # resetting some parameters for roots
    for rp in root:
        rp.hyphalEmergenceDensity = 4
        rp.lmbd = 0.15
        mycp.setOrganRandomParameter(rp)
    # setting up hyphal parameters

    # change geometry but only for hyphae
    mycp.changeGeometry(5, petri_dish)

    # check for percentage of colonized roots
    pCol = sum(mycp.getParameter("colonizationLength")) / sum(mycp.getParameter("length"))
    print("Initial colonization percentage: " + str(pCol*100) + "%")
    crossed_barrier = 0
    while pCol < 0.50:
        mycp.simulateColonization(dt,False)
        N+=1
        pCol = sum(mycp.getParameter("colonizationLength")) / sum(mycp.getParameter("length"))
        for organ in mycp.getOrgans(pb.hyphae):
            if organ.getParameter("subType") < 3:
                for node in organ.getNodes():
                    if node.x > -barrier_thickness/2 and node.z < opening_height-barrier_height and node.y < opening_length/2 and node.y > -opening_length/2:
                        crossed_barrier += 1
        if  not (crossed_barrier < 3):
            break;


    # vp.plot_roots_and_container(mycp,petri_dish)

    # print("Simulating hyphal growth until hyphae cross the barrier")
    # print(mycp.getSimTime(), max(mycp.getParameter("creationTime")))
    while crossed_barrier < 3:
        N+=1    
        # mycp.simulateHyphalGrowth(dt,False)
        # mycp.simulateHyphae(dt,False)
        mycp.simulate(dt,False)
        # print(mycp.getSimTime(), max(mycp.getParameter("creationTime")))
        # ctrs = np.array([hh.getParameter("creationTime") for hh in mycp.getOrgans(pb.root)])
        # cthhs = np.array([hh.getParameter("creationTime") for hh in  mycp.getOrgans(pb.hyphae)])
        # print(min(ctrs), max(ctrs), min(cthhs), max(cthhs))
        # print(mycp.getParameter("creationTime"))
        # raise Exception
        for organ in mycp.getOrgans(pb.hyphae):
            if organ.getParameter("subType") < 3:
                for node in organ.getNodes():
                    if node.x > -barrier_thickness/2 and node.z < opening_height-barrier_height and node.y < opening_length/2 and node.y > -opening_length/2:
                        crossed_barrier += 1

    

    # inactivating those organs that are in the root part of the compartment
    print("Inactivating those hyphae that are in the root part of the petri dish")
    mycp.turnOffSidePetriDish(-barrier_thickness/2,opening_height-barrier_height,  opening_length/2, -opening_length/2)

    # look at system to see how active
    # vp.plot_roots(mycp,"active")   
    crossed_time = mycp.getSimTime()
    # print(crossed_time)
    hours_hyphae = 60 ### HIER VERÄNDERUNG DAUER HYPHEN SIMULATION
    tip_densities = list()
    print(crossed_time,mycp.getSimTime(), max(mycp.getParameter("creationTime")))
    
    for i in range(0, hours_hyphae):
        print("Simulating hyphal growth step " + str(i+1) + " of " + str(hours_hyphae))
        # mycp.simulateHyphae(dt,False)
        mycp.simulate(dt,False)
        mycp.turnOffSidePetriDish(-barrier_thickness/2,opening_height-barrier_height,  opening_length/2, -opening_length/2)
        ana = amf.getMycSegmentAnalyser(mycp)
        if animation:
            ana.crop(small_hyphae_dish)
            ana.write("animation/" + filename + "_hoursACB_" +str(i+1) + ".vtp", ["radius", "subType", "creationTime", "organType", "colonization", "colonizationTime", "anastomosis"])
            
        # print(crossed_time, mycp.getSimTime(), max(mycp.getParameter("creationTime")))
        # raise Exception
    endsim = time.perf_counter()
    print("Time for simulation: ", endsim-start)

    # raise Exception
    # ana = amf.getMycSegmentAnalyser(mycp)
    # vp.plot_roots(mycp,"subType")

    if not animation:
        ana.write(filename + str(N+i) + ".vtp", ["radius", "subType", "creationTime", "organType", "colonization", "colonizationTime", "anastomosis","nodeTips"])

    times = np.linspace(crossed_time, max(mycp.getParameter("creationTime"))+0.01, 100)
    print('mycp.getSimTime()',mycp.getSimTime(),'max(mycp.getParameter("creationTime"))',max(mycp.getParameter("creationTime")))
    
    # times = np.linspace(0, 28, 100)

    tip_densities = amf.getParaDistperRing("nodeTips", times, ana, rings)
    times = times - crossed_time
    lengthsSubtype = amf.getParameterOverTime("length", times, ana, np.array([1,2,3]))
    return tip_densities, times, lengthsSubtype

diameter = 9.4
radius = diameter / 2

nRings = 25
location = np.array([radius*np.sqrt(i/nRings) for i in range(1, nRings+1)])
simulations = []

allsims = time.perf_counter()
for i in good_seeds:
    tip_densities, times, lengthsSubType = makesimulation(i)
    simulations.append({
    "tip_dens": tip_densities,
    "times": times[1:],
    "lengths": lengthsSubType
})
allsims_end = time.perf_counter()
print("Time for all simulations: ", allsims_end-allsims)

