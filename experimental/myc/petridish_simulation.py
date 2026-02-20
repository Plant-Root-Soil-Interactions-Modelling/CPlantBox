import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import numpy as np

mycp = pb.MycorrhizalPlant()
path = "tomatoparameters/"
name = "TomatoJohanna_WildType"

animation = True

mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)

## Setting parameters for hyphae and roots
hyphae_parameter = pb.HyphaeRandomParameter(mycp)
hyphae_parameter.subType = 1
hyphae_parameter.dx = 0.005
hyphae_parameter.a = 0.01
hyphae_parameter.b = 2.
hyphae_parameter.distTH = 0.01  # distance for anastomosis 
mycp.setOrganRandomParameter(hyphae_parameter)
# print(hyphae_parameter)

root = mycp.getOrganRandomParameter(pb.root)
for rp in root:
    rp.hyphalEmergenceDensity = 10
    rp.highresolution = 1.
    rp.dx = 0.01
    mycp.setOrganRandomParameter(rp)


## Setting up petri dish

petri_dish = pb.SDF_PlantContainer(0.94,0.94,0.1,False)
# helper_dish = pb.SDF_PlantContainer(0.94,0.94,0.02,True)
# moved_helper_dish = pb.SDF_RotateTranslate(helper_dish, 0, 0, pb.Vector3d(0.47, 0, 0))
# half_dish = pb.SDF_Difference(petri_dish, moved_helper_dish)

# vp.plot_container(moved_helper_dish)
# vp.write_container(petri_dish, "results/petri_dish.vtp")

mycp.setGeometry(petri_dish)
mycp.initialize(True)

simtime = 8
fps = 24
anim_time = simtime
N = fps * anim_time
dt = simtime / N
filename = "petri_dish_" + str(simtime)

print('First Step just for root to grow')
# infected_nodes = [0]
mycp.simulate(0.5, True)
### TODO make this better something wrong but produces some infection
print("just a check")
# for i in range(1,1000):
#     mycp.simulatePrimaryInfection(dt, True)
#     mycp.simulateSecondaryInfection(dt, True)
#     infected_nodes = mycp.getNodeInfections(2)
#     mycp.simulateHyphalGrowth(dt, True)
#     # print(infected_nodes)
# mycp.simulateHyphalGrowth(10, True)
# print(infected_nodes)
print('Rest of steps just hyphae grow')
for i in range(2, N+1):
    print('step',i, '/',N) 
    mycp.simulatePrimaryInfection(dt, True)
    mycp.simulateSecondaryInfection(dt, True)
    mycp.simulateHyphalGrowth(dt, True)     
    mycp.simulateHyphae(dt, True)

# for i in range(1, N+1):
#     print('step',i, '/',N)        
#     mycp.simulate(dt, True)
    

    
# vp.plot_plant(mycp, "organType")
vp.plot_roots_and_container(mycp,petri_dish)
# simtime = 50
# fps = 1
# anim_time = simtime
# N = fps * anim_time
# dt = simtime / N

# filename = "anastomosis_" + str(simtime)


# for i in range(1, N+1):
#     print('step',i, '/',N)
#     # print(hti)        
#     mycp.simulate(dt, False)
#     if animation:
#         ana = pb.SegmentAnalyser(mycp)
#         ana.addData("AnastomosisPoints", mycp.getAnastomosisPoints(5))
#         ana.addData("infection", mycp.getNodeInfections(2))
#         ana.write("results/" + filename + "_"+ str(i).zfill(4) + ".vtp", ["radius", "subType", "creationTime","organType","infection","AnastomosisPoints","hyphalTreeIndex"])

# # ana = pb.SegmentAnalyser(mycp)
# # ana.addData("AnastomosisPoints", mycp.getAnastomosisPoints(5))
# # ana.write("results/" + filename + ".vtp", ["radius", "subType", "creationTime","organType","infection","AnastomosisPoints","hyphalTreeIndex"])
