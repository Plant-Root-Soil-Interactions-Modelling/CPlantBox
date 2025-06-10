import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

mycp = pb.MycorrhizalPlant()
path = "../../modelparameter/structural/rootsystem/"
names = ["Anagallis_femina_Leitner_2010", "Heliantus_Pagès_2013"]

def read_initlization(name, mycp, path="../../modelparameter/structural/rootsystem/", local=True):
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
    dispersed = True
    animation = True
    infradius = 1
    for i in range(0, len(root)):
        root[i].infradius = infradius
        if root[i].infradius != 0:
            dispersed = False
            root[i].f_inf = pb.SoilLookUpSDF(infbox, 0.99, 0.0, 0.1)
        root[i].dx = 0.05
        # root[i].realize()
        # print(root[i].subType)

    mycp.initialize(True)

def simulation(simtime, frames, fps, mycp,name,local):
    N = fps * simtime
    dt = simtime / N
    filename = "RhizosphereVisual" + "_" + name
    if local:
        filename = filename + "_local"
    else:
        filename = filename + "_dispersed"
    
    for i in range(0, N):
        mycp.simulate(dt, False)
        # mycp.simulateHyphalGrowth(dt)
        if i % frames == 0:
            print("Frame " + str(i) + " of " + str(N))
            ana = pb.SegmentAnalyser(mycp)
            ana.addData("infection", mycp.getNodeInfections(2))
            ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
            ana.write(filename + "{:04d}".format(i) + ".vtp", ["radius", "subType", "creationTime", "length", "infection", "infectionTime"])
            print("Frame " + str(i) + " of " + str(N))

