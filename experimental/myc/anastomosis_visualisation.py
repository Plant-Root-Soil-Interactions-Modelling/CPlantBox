import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import math
import visualisation.vtk_plot as vp
import numpy as np

mycp = pb.MycorrhizalPlant()
path = "tomatoparameters/"
name = "TomatoJohanna_WildType"

mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)

hyphae_parameter = pb.HyphaeRandomParameter(mycp)
hyphae_parameter.subType = 1
hyphae_parameter.dx = 0.05
hyphae_parameter.distTH = 0.1  # distance for anastomosis 
mycp.setOrganRandomParameter(hyphae_parameter)
# print(hyphae_parameter)

root = mycp.getOrganRandomParameter(pb.root)
for rp in root:
    rp.hyphalEmergenceDensity = 1
    rp.highresolution = 0
    rp.dx = 0.2


mycp.initialize(True)
# print(mycp.toString())
# mycp.writeParameters(name + "_parameters.xml", 'plant', True)

simtime = 50
fps = 2
anim_time = simtime
N = fps * anim_time
dt = simtime / N

filename = "anastomosis_" + str(simtime)


for i in range(0, N):
    print('step',i, '/',N)
    # print(hti)        
    mycp.simulate(dt, False)

ana = pb.SegmentAnalyser(mycp)
ana.addData("AnastomosisPoints", mycp.getAnastomosisPoints(5))
