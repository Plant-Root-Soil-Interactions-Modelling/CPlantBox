import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

mycp = pb.MycorrhizalPlant()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
# name = "Heliantus_Pagès_2013"
mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)

hyphae_parameter = pb.HyphaeRandomParameter(mycp)
hyphae_parameter.subType = 1
hyphae_parameter.a = 0.01
hyphae_parameter.v = 1
mycp.setOrganRandomParameter(hyphae_parameter)
# print(hyphae_parameter)

root = mycp.getOrganRandomParameter(pb.root)
for rp in root:
    rp.hyphalEmergenceDensity = 1;

infbox = pb.SDF_PlantBox(10, 10, 5)
infbox = pb.SDF_RotateTranslate(infbox, 0, 0, pb.Vector3d(0, 0, -10))
local = True
animation = True
for i in range(0, len(root)):
    if local:
        root[i].f_inf = pb.SoilLookUpSDF(infbox, 0.99, 0.0, 0.1)
    root[i].dx = 0.05

mycp.initialize(True)


simtime = 10
fps = 30
anim_time = 10
N = fps * anim_time
dt = simtime / N

filename = "infection"
if animation:
    filename = "animation"
if dispersed:
    filename = filename + "_dispersed"
else:
    filename = filename + "_local"

if animation:
    for i in range(0, N):
        mycp.simulate(dt, False)
        # mycp.simulateHyphalGrowth(dt)
        ana = pb.SegmentAnalyser(mycp)
        ana.addData("infection", mycp.getNodeInfections(2))
        ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
        ana.write(filename + "{:04d}".format(i) + ".vtp", ["radius", "subType", "creationTime", "length", "infection", "infectionTime"])
        print("Frame " + str(i) + " of " + str(N))

else:
    for i in range(0, N):
        mycp.simulate(dt, False)
        # roost = mycp.getOrganRandomParameter(pb.root)
        # for i in range(0, len(roost)):
        #     print(root[i].subType)
        # mycp.simulateHyphalGrowth(dt)

    # mycp.simulate(simtime, False)
    # print("sim time", mycp.getSimTime())

    # mycp.simulateHyphalGrowth(simtime)
    # hyphae = mycp.getOrgans(5)
    # print("number of hyphae", len(hyphae))
    # print("type", type(hyphae))
    # for h in hyphae:
    #     print(h.getParameter("age"))
    

    ana = pb.SegmentAnalyser(mycp)
    ana.addData("infection", mycp.getNodeInfections(2))
    ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
    pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "length", "infection", "infectionTime","organType"])
    vp.plot_roots(ana, "infection")
    # vp.plot_roots(ana, "infectionTime")
    # vp.plot_plant(mycp, "organType")
    ana.write(filename + ".vtp", ["radius", "subType", "creationTime", "length", "infection", "infectionTime","organType"])

