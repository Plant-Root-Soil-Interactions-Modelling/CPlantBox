import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import numpy as np

mycp = pb.MycorrhizalPlant()
path = "../tomatoparameters/"
name = "TomatoJohanna_WildType"

plant = pb.Plant()


animation = True

mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)
plant.readParameters(path + name + ".xml", fromFile = True, verbose = True)

## Setting parameters for hyphae and roots
hyphae_parameter = pb.HyphaeRandomParameter(mycp)
hyphae_parameter.subType = 1
hyphae_parameter.dx = 0.001
hyphae_parameter.a = 0.01
hyphae_parameter.b = 2.
hyphae_parameter.distTH = 0.01  # distance for anastomosis 
mycp.setOrganRandomParameter(hyphae_parameter)
# print(hyphae_parameter)

root = mycp.getOrganRandomParameter(pb.root)
for rp in root:
    rp.hyphalEmergenceDensity = 0
    rp.highresolution = 0
    rp.dx = 0.01
    rp.a = 0.01
    mycp.setOrganRandomParameter(rp)

root = plant.getOrganRandomParameter(pb.root)
for rp in root:
    rp.dx = 0.01
    rp.a = 0.01
    plant.setOrganRandomParameter(rp)

# make sure to set the seed position to 0 because of the petri dish
seed_parameter = pb.SeedRandomParameter(mycp)
seed_parameter.seedPos.z = 0
seed_parameter.seedPos.x = 0
seed_parameter.seedPos.y = 0
mycp.setOrganRandomParameter(seed_parameter)

seed_parameter = pb.SeedRandomParameter(plant)
seed_parameter.seedPos.z = 0
seed_parameter.seedPos.x = 0
seed_parameter.seedPos.y = 0
plant.setOrganRandomParameter(seed_parameter)

## Setting up petri dish

petri_dish = pb.SDF_PlantContainer(0.94,0.94,0.1,False)
helper_dish = pb.SDF_PlantContainer(0.94,0.94,0.1,True)
moved_helper_dish = pb.SDF_RotateTranslate(helper_dish, 0, 0, pb.Vector3d(-0.94, 0, 0))
half_dish = pb.SDF_Difference(petri_dish, moved_helper_dish)

# vp.plot_container(half_dish)
# vp.write_container(petri_dish, "results/petri_dish.vtp")

mycp.setGeometry(half_dish)
plant.setGeometry(half_dish)
mycp.initialize(True)
plant.initialize(True)

simtime = 15
fps = 24
anim_time = simtime
N = fps * anim_time
dt = simtime / N

filename = "half_dish_" + str(simtime)


for i in range(0, N):
    if i % 6 == 0:
        print("Step " + str(i) + " of " + str(N))
    mycp.simulate(dt,True)

plant.simulate(simtime,True)
vp.plot_roots_and_container(mycp,half_dish)
vp.plot_roots_and_container(plant,half_dish)
mycp