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
fpd = 24
N = simtime * fpd
dt = simtime / N
time = np.linspace(0, simtime, N)

primInfL =[] # legnth of primary infection
secInfL = [] # legnth of secondary infection
lenghtL = [] # legnth of the root system
nonmycL = [] # legnth of the non-mycorrhizal part
ratioL = [] # ratio of primary infection to non-mycorrhizal part from each time step

mycp.simulate(dt, False)
primInfL.append(sum(mycp.getParameter("primaryInfection")))
secInfL.append(sum(mycp.getParameter("secondaryInfection")))
lenghtL.append(sum(mycp.getParameter("length")))
nonmycL.append(lenghtL[-1]-primInfL[-1]-secInfL[-1])
ratioL.append(abs(primInfL[-1])/(nonmycL[-1]))

for t in range(1,N):
    mycp.simulate(dt, False)
    primInfL.append(sum(mycp.getParameter("primaryInfection")))
    secInfL.append(sum(mycp.getParameter("secondaryInfection")))
    lenghtL.append(sum(mycp.getParameter("length")))
    nonmycL.append(lenghtL[-2]-primInfL[-2]-secInfL[-2])
    ratioL.append(abs(primInfL[-1]-primInfL[-2])/(nonmycL[-1]))

ratio = True
if ratio:
    plt.plot(time, np.asarray(ratioL), label="Primary Infection Ratio")
    plt.title("Infection Ratio over time")
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Ratio")
    plt.show()
else:
    plt.plot(time, np.asarray(primInfL), label="Primary Infection")
    plt.plot(time, np.asarray(secInfL), label="Secondary Infection")
    plt.plot(time, np.asarray(lenghtL), label="Length")
    plt.legend()
    plt.title("Infection over time")
    plt.xlabel("Time")
    plt.ylabel("[cm]")
    plt.show()