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

## Setting up petri dish
diameter = 9.4
radius = diameter / 2
height = 1.6

# petri dish has a radius of 9.4 cm and a height of 1.6 cm
petri_dish = pb.SDF_PlantContainer(radius,radius,height,False)
# the helper dish is used to cut the petri dish in half, it has the same radius and height as the petri dish but is rotated and translated to cut the petri dish in half
helper_dish = pb.SDF_PlantContainer(radius,radius,height,True)

moved_helper_dish = pb.SDF_RotateTranslate(helper_dish, 0, 0, pb.Vector3d(-radius, 0, 0))
half_dish = pb.SDF_Difference(petri_dish, moved_helper_dish)
# vp.plot_container(half_dish)

# introduce parameters for barrier and opening
barrier_thickness = 0.16
barrier_height = height
opening_length = 5.0
opening_height = 0.2

helper_staff = pb.SDF_PlantBox(barrier_thickness,barrier_height,diameter)
helper_staff2 = pb.SDF_PlantBox(barrier_thickness,opening_height,opening_length)

# have to  move the helper staff for the right position of the opening and barrier, the midpoint is the distance from the center of the helper staff to the center of the petri dish, the bottompoint is the distance from the center of the helper staff to the center of the petri dish in the y direction, and the helper staff is moved to the position of the opening and barrier
midpoint = opening_length / 2 - diameter / 2
bottompoint = opening_height / 2 - barrier_height / 2
helper_staff2 = pb.SDF_RotateTranslate(helper_staff2, 0,0,pb.Vector3d(0, bottompoint, midpoint))
helper_dish2 = pb.SDF_Difference(helper_staff, helper_staff2)
moved_helper_dish2 = pb.SDF_RotateTranslate(helper_dish2, 90, pb.SDF_Axis.xaxis , pb.Vector3d(0, -radius, -height/2))
petri_dish = pb.SDF_Difference(petri_dish, moved_helper_dish2)

root = mycp.getOrganRandomParameter(pb.root)
for rp in root:
    rp.hyphalEmergenceDensity = 0
    rp.highresolution = 0
    rp.dx = 0.1
    # rp.a = 0.01
    rp.tropismT = 2
    mycp.setOrganRandomParameter(rp)


# make sure to set the seed position to 0 because of the petri dish
seed_parameter = pb.SeedRandomParameter(mycp)
seed_parameter.seedPos.z = -height / 2 # seed is positioned in the middle of the petri dish in the z direction
seed_parameter.seedPos.x = 2.5
seed_parameter.seedPos.y = 0
mycp.setOrganRandomParameter(seed_parameter)
mycp.setGeometry(half_dish)
mycp.initialize(True)

simtime = 10
fps = 24
anim_time = simtime
N = fps * anim_time
dt = simtime / N

filename = "half_dish_" + str(simtime)

start = time.perf_counter()

for i in range(0, N):
    if i % 6 == 0:
        print("Step " + str(i) + " of " + str(N))
    mycp.simulate(dt,True)

middle = time.perf_counter()
for rp in root:
    rp.hyphalEmergenceDensity = 1
    rp.highresolution = 1.
    mycp.setOrganRandomParameter(rp)

# Setting parameters for hyphae and roots
hyphae_parameter = mycp.getOrganRandomParameter(pb.hyphae)
for hp in hyphae_parameter:    
    hp.a = 0.01
    hp.ln = 0.1
    hp.b_prob = 0.
    # hp.lmax = 0.5 
    hp.lb = 0.01
    # hp.dx = 0.001
    hp.v = 0.13
    hp.tropismS = 0.1
    hp.distTH = 0.01   # distance for anastomosis
    mycp.setOrganRandomParameter(hp)

mycp.changeGeometry(5, petri_dish)
# need something here to redo the tropism function just for the hyphae
dt_hyphalgrowth = 1/24/60
print("Maximal distance a hypha can grow in one time step: " + str(dt_hyphalgrowth * hyphae_parameter[0].v) + " cm")
print("Anastomosis distance: " + str(hyphae_parameter[0].distTH) + " cm")
mycp.simulateHyphalGrowth(1,False)
middelsimHyphae = time.perf_counter()
for i in range(0, 9):
    mycp.simulateHyphae(dt_hyphalgrowth,False)
end = time.perf_counter()

print(f"Time for root growth: {middle - start:.2f} seconds")
print(f"Time for hyphal creation: {end - middle:.2f} seconds")
print(f"Time for hyphal simulation: {end - middelsimHyphae:.2f} seconds")
# second dish just to see hyphal growth
# half_dish2 = pb.SDF_TranslateRotate(half_dish, 0, 0, pb.Vector3d(1.5, 0, 0))
vp.plot_roots_and_container(mycp,petri_dish)
ana = pb.SegmentAnalyser(mycp)
ana.addData("infection", mycp.getNodeInfections(2))
ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
ana.addData("anastomosis", mycp.getAnastomosisPoints(5))
ana.write(filename + ".vtp", ["radius", "subType", "creationTime", "organType", "infection", "infectionTime", "anastomosis"])