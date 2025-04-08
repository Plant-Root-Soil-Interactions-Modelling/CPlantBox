import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import matplotlib.pyplot as plt
import numpy as np


mycp = pb.MycorrhizalPlant()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
# name = "Heliantus_Pagès_2013"
mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)

# plant = pb.Organism()  # store organism (not owned by Organ, or OrganRandomParameter)
# p0 = pb.MycorrhizalRootRandomParameter(plant)
# p0.name, p0.subType, p0.la, p0.lb, p0.lmax, p0.ln, p0.r, p0.dx, p0.infradius = "taproot", 1, 10., 1., 100., 1., 1.5, 0.5, 0.0
# p0.successor = [[2]]
# p0.successorP = [[1.]]
# p1 = pb.MycorrhizalRootRandomParameter(plant)
# p1.name, p1.subType, p1.lmax, p1.r, p1.dx, p1.infradius = "lateral", 2, 25., 2., 0.1, 0.0
# p0, p1 = p0, p1  # needed at later point
# plant.setOrganRandomParameter(p0)  # the organism manages the type parameters and takes ownership
# plant.setOrganRandomParameter(p1)
# srp = pb.SeedRandomParameter(plant)
# plant.setOrganRandomParameter(srp)

# param0 = p0.realize()  # set up root by hand (without a root system)
# param0.la, param0.lb = 0, 0  # its important parent has zero length, otherwise creation times are messed up
# parentroot = pb.MycorrhizalRoot(1, param0, True, True, 0., 0., pb.Vector3d(0, 0, -1), 0,False, 0)  # takes ownership of param0
# parentroot.setOrganism(plant)
# parentroot.addNode(pb.Vector3d(0, 0, -3), 0)  # there is no nullptr in Python

# parentroot = parentroot  # store parent (not owned by child Organ)
# mycroot = pb.MycorrhizalRoot(plant, p0.subType,  0, parentroot, 0)
# mycroot.setOrganism(plant)

hyphae_parameter = pb.HyphaeRandomParameter(mycp)
hyphae_parameter.subType = 1
hyphae_parameter.a = 0.01
hyphae_parameter.v = 1
mycp.setOrganRandomParameter(hyphae_parameter)
# print(hyphae_parameter)

root = mycp.getOrganRandomParameter(pb.root)
for rp in root:
    rp.hyphalEmergenceDensity = 1;

infbox = pb.SDF_PlantBox(3, 3, 3)
infbox = pb.SDF_RotateTranslate(infbox, 0, 0, pb.Vector3d(0, 0, -10))
dispersed = True
animation = False
infradius = 0
for i in range(0, len(root)):
    root[i].infradius = infradius
    if root[i].infradius != 0:
        dispersed = False
    root[i].dx = 0.05
    root[i].f_inf = pb.SoilLookUpSDF(infbox, 1, 0.0, 0.1)

mycp.initialize(True)


simtime = 10

primInfL =[]
secInfL = []
lenghtL = []
ratioL = []


for t in range(0,simtime):
    mycp.simulate(1, False)
    # print(mycp.getParameter("primaryInfection"))
    primInfL.append(sum(mycp.getParameter("primaryInfection")))
    secInfL.append(sum(mycp.getParameter("secondaryInfection")))
    lenghtL.append(sum(mycp.getParameter("length")))
    ratioL.append(sum(mycp.getParameter("primaryInfection"))/sum(mycp.getParameter("length")))
ratio = True
if ratio:
    plt.plot(np.asarray(ratioL), label="Primary Infection Ratio")
    plt.title("Infection Ratio over time")
    plt.xlabel("Time")
    plt.ylabel("[cm]")
    plt.show()
else:
    plt.plot(np.asarray(primInfL), label="Primary Infection")
    plt.plot(np.asarray(secInfL), label="Secondary Infection")
    plt.plot(np.asarray(lenghtL), label="Length")
    plt.title("Infection over time")
    plt.xlabel("Time")
    plt.ylabel("[cm]")
    plt.show()