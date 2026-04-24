import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import numpy as np

import time

mycp = pb.MycorrhizalPlant()
path = "../tomatoparameters/"
name = "TomatoJohanna_WildTypeTwoHyphaeTypes"

animation = True
mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)

### initial root parameters
root = mycp.getOrganRandomParameter(pb.root)
for rp in root:
    rp.hyphalEmergenceDensity = 0
    rp.highresolution = 0.
    rp.dx = 0.1
    rp.maxAge = 100 # maximal infection age
    rp.hyphalDelay = 5.0
    # rp.a = 0.01
    rp.tropismT = 1
    mycp.setOrganRandomParameter(rp)


## Setting up petri dish
diameter = 9.4
radius = diameter / 2
height = 1.6

# introduce parameters for barrier and opening
barrier_thickness = 0.16
barrier_height = height
opening_length = 5.0
opening_height = 0.2

# petri dish has a radius of 9.4 cm and a height of 1.6 cm
petri_dish = pb.SDF_PlantContainer(radius,radius,height,False)
# the helper dish is used to cut the petri dish in half, it has the same radius and height as the petri dish but is rotated and translated to cut the petri dish in half
helper_dish = pb.SDF_PlantContainer(radius,radius,height,True)
# moving the helper dish such that it halves the petri dish and removes a bit more to restrict roots to the correct side of barrier
moved_helper_dish = pb.SDF_RotateTranslate(helper_dish, 0, 0, pb.Vector3d(-(radius+barrier_thickness+root[0].a), 0, 0))
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

moved_helper_dish_hyphae = pb.SDF_RotateTranslate(helper_dish, 0, 0, pb.Vector3d(-(radius+barrier_thickness), 0, 0))
hyphae_petri_dish = pb.SDF_Difference(petri_dish, moved_helper_dish_hyphae)

hyphae_dishes = []
for i in range(1, 10):
    small_dish = pb.SDF_PlantContainer(radius*i/10,radius*i/10,height,False)
    small_hyphae_dish = pb.SDF_Intersection(small_dish, moved_helper_dish_hyphae)
    hyphae_dishes.append(small_hyphae_dish)

# make sure to set the seed position to 0 because of the petri dish
seed_parameter = pb.SeedRandomParameter(mycp)
seed_parameter.seedPos.z = -height / 6 # seed is positioned in the middle of the petri dish in the z direction
seed_parameter.seedPos.x = - 1.0
seed_parameter.seedPos.y = 0
mycp.setOrganRandomParameter(seed_parameter)

mycp.setGeometry(half_dish)
mycp.initialize(True)

# set up simulation times etc.
simtime = 20
fps = 24
anim_time = simtime
N = fps * anim_time
dt = simtime / N

# Start simulation
start = time.perf_counter()
for i in range(0, N):
    if i % 6 == 0:
        print("Step " + str(i) + " of " + str(N))
    mycp.simulate(dt,True)
    if (animation):
        ana = pb.SegmentAnalyser(mycp)
        ana.addData("infection", mycp.getNodeInfections(2))
        ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
        ana.addData("anastomosis", mycp.getAnastomosisPoints(5))
        ana.write("animation/step" + str(i) + ".vtp", ["radius", "subType", "creationTime", "organType", "infection", "infectionTime", "anastomosis"])
# look at roots and container
vp.plot_roots_and_container(mycp,half_dish)


afterroots = time.perf_counter()
# resetting some parameters for roots
for rp in root:
    rp.hyphalEmergenceDensity = 2
    mycp.setOrganRandomParameter(rp)
# setting up hyphal parameters
hyphae_parameter = mycp.getOrganRandomParameter(pb.hyphae)
for hp in hyphae_parameter:    
    hp.a = 0.01
    hp.ln = 0.1
    hp.b_prob = 0.
    hp.lb = 0.01
    hp.dx = 0.01
    hp.v = 0.13
    hp.tropismS = 1.0
    hp.distTH = 0.01   # distance for anastomosis
    mycp.setOrganRandomParameter(hp)

# change geometry but only for hyphae
mycp.changeGeometry(5, petri_dish)

# check for percentage of colonized roots
pCol = sum(mycp.getParameter("infectionLength")) / sum(mycp.getParameter("length"))
print("Initial infection percentage: " + str(pCol*100) + "%")
days = 0
while pCol < 0.60:
    mycp.simulateInfection(1,False)
    days += 1
    N+=24
    pCol = sum(mycp.getParameter("infectionLength")) / sum(mycp.getParameter("length"))
    print("Infection percentage: " + str(pCol*100) + "% after " + str(days) + " days.")
print("Infection reached 60% after " + str(days) + " days.")

# simulating hyphal growth


mycp.simulateHyphalGrowth(0.5,False)
vp.plot_roots_and_container(mycp,petri_dish)
# print("Simulating hyphal growth until hyphae cross the barrier")
# crossed_barrier = False

# while not crossed_barrier:
#     N+=12    
#     mycp.simulateHyphalGrowth(dt,False)
#     mycp.simulateHyphae(dt,False)
#     for organ in mycp.getOrgans():
#         for node in organ.getNodes():
#             if node.x > -barrier_thickness/2 and node.z < barrier_height-opening_height:
#                 crossed_barrier = True
#     if (N % 20 == 0):
#         vp.plot_roots_and_container(mycp,petri_dish)
#     ana = pb.SegmentAnalyser(mycp)
#     ana.addData("infection", mycp.getNodeInfections(2))
#     ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
#     ana.addData("anastomosis", mycp.getAnastomosisPoints(5))
#     ana.write("animation/step" + str(N+1) + ".vtp", ["radius", "subType", "creationTime", "organType", "infection", "infectionTime", "anastomosis"])


# inactivating those organs that are in the root part of the compartment
organs = mycp.getOrgans()
for organ in organs:
    stayactive = False
    for node in organ.getNodes():
        if node.x > -barrier_thickness/2 and node.z < barrier_height-opening_height:
            stayactive = True
    organ.setActive(stayactive)

# look at system to see how active
vp.plot_roots(mycp,"active")   

middelsimHyphae = time.perf_counter()

for i in range(1, 60):
    print("Simulating hyphal growth step " + str(i+1) + " of 60")
    mycp.simulateHyphae(dt,False)
    ana = pb.SegmentAnalyser(mycp)
    ana.addData("infection", mycp.getNodeInfections(2))
    ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
    ana.addData("anastomosis", mycp.getAnastomosisPoints(5))
    ana.write("animation/step" + str(N+i) + ".vtp", ["radius", "subType", "creationTime", "organType", "infection", "infectionTime", "anastomosis"])

end = time.perf_counter()
filename = "splitpetri_dish" + str(simtime)

vp.plot_roots_and_container(mycp,petri_dish)
vp.write_container(petri_dish, "petri_dish.vtp")