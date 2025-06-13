import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

mycp = pb.MycorrhizalPlant()
path = "../../modelparameter/structural/rootsystem/"
names = ["Anagallis_femina_Leitner_2010", "Heliantus_Pagès_2013","maize","Glycine_max","Faba_synMRI","Brassica_napus_a_Leitner_2010","Pisum_sativum_a_Pagès_2014"]

# --- Begin code previously in read_initlization ---
name = names[3]
local = True

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

infbox = pb.SDF_PlantBox(10, 10, 10)
infbox = pb.SDF_RotateTranslate(infbox, 0, 0, pb.Vector3d(0, 0, -20))

for i in range(0, len(root)):
    if local:
        root[i].f_inf = pb.SoilLookUpSDF(infbox, 0.99, 0.0, 0.1)
    root[i].dx = 0.05
    # root[i].realize()
    # print(root[i].subType)

mycp.initialize(True)
# --- End code previously in read_initlization ---

# --- Begin code previously in simulation ---
simtime = 20
fps = 15
N = fps * simtime
dt = simtime / N
filename = "RhizosphereVisual" + "_" + name
if local:
    filename = filename + "_local_"
else:
    filename = filename + "_dispersed_"

for i in range(0, N):
    mycp.simulate(dt, False)
    # mycp.simulateHyphalGrowth(dt)
    if i % (fps*2) == 0:
        print("Frame " + str(i) + " of " + str(N))
        ana = pb.SegmentAnalyser(mycp)
        ana.addData("infection", mycp.getNodeInfections(2))
        ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
        ana.write(filename + "{:04d}".format(i) + ".vtp", ["radius", "subType", "creationTime", "length", "infection", "infectionTime"])
        print("Frame " + str(i) + " of " + str(N))
# --- End code previously in simulation ---

