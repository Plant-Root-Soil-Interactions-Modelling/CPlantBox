import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import numpy as np

# set plant
mycp = pb.MycorrhizalPlant(100)
path = "tomatoparameters/"
name = "TwoHyphaePlusBAS"
# set geomtry of pot check Henri´s pot sizes
depth = 30
radius = 11.7
pot = pb.SDF_PlantContainer(radius, radius, depth, False)

mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)

root = mycp.getOrganRandomParameter(pb.root)
for rp in root:
    rp.hyphalEmergenceDensity = 2
    rp.highresolution = 0
    rp.dx = 0.1
    # rp.a = 0.01
    rp.tropismT = 2
    mycp.setOrganRandomParameter(rp)

mycp.setGeometry(pot)
animation = True

seed_parameter = pb.SeedRandomParameter(mycp)
seed_parameter.seedPos.z = -3 # seed is positioned in the middle of the petri dish in the z direction
seed_parameter.seedPos.x = 0
seed_parameter.seedPos.y = 0
mycp.setOrganRandomParameter(seed_parameter)

simtime = 25
fps = 24
anim_time = simtime
N = fps * anim_time
dt = simtime / N

filename = "potsim_" + str(simtime/24)
mycp.initialize()
for i in range(0, N):
    if i % 6 == 0:
        print("Step " + str(i) + " of " + str(N))
    mycp.simulate(dt,False)
    if (animation):
        ana = pb.SegmentAnalyser(mycp)
        ana.addData("infection", mycp.getNodeInfections(2))
        ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
        ana.addData("anastomosis", mycp.getAnastomosisPoints(5))
        ana.write("animation/" + filename + "_" + str(i) + ".vtp", ["radius", "subType", "creationTime", "organType", "infection", "infectionTime", "anastomosis"])
# set parameters for roots and hyphae
vp.write_container(pot, "pot.vtp")
# run simulation and save results

# extract soil cores

# analsyse soil cores for hyphal density