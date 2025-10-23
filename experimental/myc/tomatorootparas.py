import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import math
import visualisation.vtk_plot as vp
import numpy as np

expWSLength = 9229 # cm
expWLength = 7303 # cm

expRMCSLength = 14301 # cm
expRMCLength = 10639 # cm

expWSDia = 0.33
expWDia = 0.34

expRMCSDia = 0.32
expRMCDia = 0.27


mycp = pb.MycorrhizalPlant()
# path = "../../modelparameter/structural/rootsystem/"
path = ""
name = "TomatoJohanna"

mycp.readParameters(path +name + ".xml", fromFile = True, verbose = True)

hyphae_parameter = pb.HyphaeRandomParameter(mycp)
hyphae_parameter.subType = 1
hyphae_parameter.a = 0.01
hyphae_parameter.dx = 0.01
mycp.setOrganRandomParameter(hyphae_parameter)


root = mycp.getOrganRandomParameter(pb.root)
for rp in root:
    rp.hyphalEmergenceDensity = 0
    rp.highresolution = 0
    rp.dx = 0.2

mycp.initialize(True)
totalduration = 108
potduration = 73
simtime = totalduration
fpd = 1
N = simtime * fpd
dt = simtime / N

for i in range(1,N):
        print('step',i, '/',N)
        mycp.simulate(dt, False)
        # mycp.simulateHyphalGrowth(dt)
    
print('done')
vp.plot_plant(mycp, "organType")
ana = pb.SegmentAnalyser(mycp)
rad = ana.getParameter("radius")
rad = sum(rad)/len(rad)
rl = ana.getSummed("length")
print("Root length total (cm):", rl)
print("Average root diameter (cm):", rad)
# ana.addData("infection", mycp.getNodeInfections(2))
# ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
# pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "length", "infection", "infectionTime","organType"])
vp.plot_roots(ana, "radius")
# # vp.plot_roots(ana, "infectionTime")
# # vp.plot_plant(mycp, "organType")
# ana.write(name + ".vtp", ["radius", "subType", "creationTime","organType"])# "infection", "infectionTime",