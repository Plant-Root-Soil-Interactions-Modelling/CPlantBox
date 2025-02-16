import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp


mycp = pb.MycorrhizalPlant(0)
p = pb.Plant(0)
path = "../../modelparameter/structural/rootsystem/"
name = "Heliantus_Pagès_2013"
mycp.readParameters(path + name + ".xml",fromFile = True, verbose = True)

mycp.initialize(True)
# print(mycp.getOrganRandomParameter(2))
mycp.simulate(30, True)
ana = pb.SegmentAnalyser(mycp)
ana.addData("infection",mycp.getNodeInfections(2))
ana.addData("infectionTime", mycp.getNodeIT(2))
pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "length", "infection", "infectionTime"])
vp.plot_roots(pd, "infectionTime") 