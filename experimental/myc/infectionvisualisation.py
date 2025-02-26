import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp


mycp = pb.MycorrhizalPlant()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
# name = "Heliantus_Pagès_2013"
mycp.readParameters(path + name + ".xml",fromFile = True, verbose = True)

mycp.initialize(True)
# print(mycp.getOrganRandomParameter(2))

simtime = 30



mycp.simulate(simtime, True)
ana = pb.SegmentAnalyser(mycp)
ana.addData("infection",mycp.getNodeInfections(2))
ana.addData("infectionTime", mycp.getNodeIT(2))
pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "length", "infection", "infectionTime"])
vp.plot_roots(pd, "infection") 