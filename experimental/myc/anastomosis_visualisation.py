import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
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
    mycp.setOrganRandomParameter(rp)

mycp.initialize(True)
# print(mycp.toString())
# mycp.writeParameters(name + "_parameters.xml", 'plant', True)

simtime = 50
fps = 1
anim_time = simtime
N = fps * anim_time
dt = simtime / N

filename = "anastomosis_" + str(simtime)


for i in range(1, N+1):
    print('step',i, '/',N)
    # print(hti)        
    mycp.simulate(dt, False)

print(mycp.getAnastomosisPoints(5)[0])
points = pb.SegmentPoints(mycp.getAnastomosisPoints(5))

points.write("anastomosis_points.vtp", True)
# ana = pb.SegmentAnalyser(mycp)
# ana.addData("AnastomosisPoints", mycp.getAnastomosisPoints(5)) TODO Does not work bc addData does not support 3D Vectors
# ana.write(filename + ".vtp", True)

