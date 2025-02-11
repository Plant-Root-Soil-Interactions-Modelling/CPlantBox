import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp


mycp = pb.MycorrhizalPlant()
path = "../../modelparameter/structural/rootsystem/"
name = "Heliantus_Pagès_2013"
mycp.readParameters(path + name + ".xml",fromFile = True, verbose = True)

mycp.initialize()
# print(mycp.getOrganRandomParameter(2))
mycp.simulate(30, True)
params = mycp.getOrganRandomParameter(2)
print(params[1])

CTs = mycp.getNodeCTs()
print(CTs)

# infs = mycp.getNodeInfections(pb.root)
# print(infs)

# ana = pb.SegmentAnalyser(mycp)
# ana.addData("Infection",mycp.getNodeInfections(2))
# ana.addData("rx", rx_b)  # node data are converted to segment data
# ana.addData("radial", radial_fluxes)
# ana.addData("axial", axial_fluxes)
# ana.addData("net", axial_i-axial_j-radial_fluxes) # np.maximum(1.e-2,)
# pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "length", "rx", "axial", "radial", "net"])
# vp.plot_roots(pd, "net") # axial, radial, rx