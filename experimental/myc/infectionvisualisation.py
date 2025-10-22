import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import math
import visualisation.vtk_plot as vp
import numpy as np

mycp = pb.MycorrhizalPlant(1)
path = "../../modelparameter/structural/rootsystem/"
name = "Glycine_max"
# name = "Heliantus_Pagès_2013"

animation = False
local = True
infbox = pb.SDF_PlantBox(4, 4, 4)
# infbox = pb.SDF_RotateTranslate(infbox, 0, 0, pb.Vector3d(0, 0, -10))
for i in range(1,5):
    infbox = pb.SDF_Union(infbox, pb.SDF_RotateTranslate(infbox, 0, 0, pb.Vector3d(0, 0, - i * 10)))

mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)

hyphae_parameter = pb.HyphaeRandomParameter(mycp)
hyphae_parameter.subType = 1
hyphae_parameter.a = 0.01
hyphae_parameter.dx = 0.01
hyphae_parameter.distTH = 0.001  # distance for anastomosis
mycp.setOrganRandomParameter(hyphae_parameter)
# print(hyphae_parameter)

root = mycp.getOrganRandomParameter(pb.root)
for rp in root:
    rp.hyphalEmergenceDensity = 2
    rp.highresolution = 1
    if local:
        rp.f_inf = pb.SoilLookUpSDF(infbox, 0.99, 0.0, 0.1)
    rp.dx = 0.2
    # rp.hyphalEmergenceDist = 0.5;  # distance between hyphal emergence points


mycp.initialize(True)
# print(mycp.toString())
# mycp.writeParameters(name + "_parameters.xml", 'plant', True)

simtime = 10
fps = 30
anim_time = 10
N = fps * anim_time
dt = simtime / N

filename = "infection_" + str(simtime)
if animation:
    filename = "animation"
if not local:
    filename = filename + "_dispersed"
else:
    filename = filename + "_local"


if animation:
    for i in range(1, N):
        mycp.simulate(dt, False)
        # mycp.simulateHyphalGrowth(dt)
        ana = pb.SegmentAnalyser(mycp)
        ana.addData("infection", mycp.getNodeInfections(2))
        ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
        ana.write(filename + "{:04d}".format(i) + ".vtp", ["radius", "subType", "creationTime", "length", "infection", "infectionTime"])
        print("Frame " + str(i) + " of " + str(N))

else:
    for i in range(0, N):
        print('step',i, '/',N)
        mycp.simulate(dt, False)
        # mycp.simulateHyphalGrowth(dt)
    
    print('done')
    
    ana = pb.SegmentAnalyser(mycp)
    ana.addData("infection", mycp.getNodeInfections(2))
    ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
    pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "length", "infection", "infectionTime","organType"])
    vp.plot_roots(ana, "infection")
    # vp.plot_roots(ana, "infectionTime")
    # vp.plot_plant(mycp, "organType")
    ana.write(filename + ".vtp", ["radius", "subType", "creationTime","organType"])# "infection", "infectionTime",

