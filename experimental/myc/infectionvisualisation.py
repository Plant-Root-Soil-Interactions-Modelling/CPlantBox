import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp


mycp = pb.MycorrhizalPlant()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
# name = "Heliantus_Pagès_2013"
mycp.readParameters(path + name + ".xml",fromFile = True, verbose = True)

root = mycp.getOrganRandomParameter(pb.root)
infbox = pb.SDF_PlantBox(3, 3, 3)
infbox = pb.SDF_RotateTranslate(infbox, 0, 0, pb.Vector3d(0, 0, -10))
dispersed = True
animation = False
infradius = 1
for i in range(0,len(root)):
    root[i].infradius = infradius
    if root[i].infradius != 0:
        dispersed = False
    root[i].dx = 0.1
    root[i].f_inf = pb.SoilLookUpSDF(infbox, 1, 0.0, 0.5)
    # print(root[i].f_inf.getValue(pb.Vector3d(0., 0., -11.)))

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
    for i in range(0,N):
        mycp.simulate(dt, False)
        ana = pb.SegmentAnalyser(mycp)
        ana.addData("infection",mycp.getNodeInfections(2))
        ana.addData("infectionTime", mycp.getNodeIT(2))
        ana.write( filename + "{:04d}".format(i)+".vtp",["radius", "subType", "creationTime", "length", "infection", "infectionTime"])
        print("Frame " + str(i) + " of " + str(N))

else:
    mycp.simulate(simtime, True)
    ana = pb.SegmentAnalyser(mycp)
    ana.addData("infection",mycp.getNodeInfections(2))
    ana.addData("infectionTime", mycp.getNodeIT(2))
    pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "length", "infection", "infectionTime"])
    vp.plot_roots(pd, "infection") 
    ana.write("nonanimation_local.vtp",["radius", "subType", "creationTime", "length", "infection", "infectionTime"])