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
infradius = 0
for i in range(0, len(root)):
    root[i].infradius = infradius
    if root[i].infradius != 0:
        dispersed = False
        print("ATTENTION: The infection radius is not 0, the infection will be local")
    root[i].dx = 0.05
    root[i].f_inf = pb.SoilLookUpSDF(infbox, 1, 0.0, 0.1)

mycp.initialize(True)


simtime = 10
fpd = 24
N = simtime * fpd
dt = simtime / N
time = np.linspace(0, simtime, N)

observed_primary = [] # length of primary infection
observed_new_primary = [] # length of new primary infection
expected_new_primary = [] # expected new primary infection
observed_secondary = [] # length of secondary infection
total_lengths = [] # length of the root system
# nonmycL = [] # length of the non-mycorrhizal part
# ratioL = [] # ratio of primary infection to non-mycorrhizal part from each time step

mycp.simulate(dt, False)
observed_primary.append(sum(mycp.getParameter("primaryInfection")))
observed_secondary.append(sum(mycp.getParameter("secondaryInfection")))
total_lengths.append(sum(mycp.getParameter("length")))
# nonmycL.append(total_lengths[-1]-observed_new_primary[-1]-observed_secondary[-1])
# ratioL.append(abs(observed_new_primary[-1])/(nonmycL[-1])) 


for t in range(1, len(time)):
    mycp.simulate(dt, False)
    observed_primary.append(sum(mycp.getParameter("primaryInfection")))
    observed_secondary.append(sum(mycp.getParameter("secondaryInfection")))
    total_lengths.append(sum(mycp.getParameter("length")))

    # Differenzen
    delta_obs = observed_primary[t] - observed_primary[t-1]
    observed_new_primary.append(delta_obs)

    # Nicht infizierte Länge zum Start dieses Zeitschritts
    L_susceptible = (total_lengths[t-1]
                     - observed_primary[t-1]
                     - observed_secondary[t-1])

    delta_expected =  mycp.getOrganRandomParameter(pb.root)[1].p * L_susceptible * dt
    expected_new_primary.append(delta_expected)

    # if t % 10 == 0:
    #     print(f"t={t}:")
    #     print(f"- Beobachtet:   {delta_obs:.2f} mm neue primäre Infektion")
    #     print(f"- Erwartet:     {delta_expected:.2f} mm")
    #     print(f"- Abweichung:   {delta_obs - delta_expected:.2f} mm")
    #     print(f"- Totale Länge: {total_lengths[t]:.2f} mm")
    #     print(f" - Nicht infizierte Länge: {L_susceptible:.2f} mm")
    #     print(f"- Primäre Infektion: {observed_primary[t]:.2f} mm")
    #     print(f"- Sekundäre Infektion: {observed_secondary[t]:.2f} mm")
    #     print("-" * 30)

ratio = True
if ratio:
    plt.figure(figsize=(10, 6))
    plt.plot(time[1:], observed_new_primary, label="Primary Infection")
    plt.plot(time[1:], expected_new_primary, label="Expected Primary Infection")
    plt.title("Expected vs Observed Primary Infection")
    plt.legend()
    plt.xlabel("Time Step")
    plt.ylabel("Primary Infection Length [cm]")
    plt.show()
else:
    plt.plot(np.asarray(observed_primary), label="Primary Infection")
    plt.plot(np.asarray(observed_secondary), label="Secondary Infection")
    plt.plot(np.asarray(total_lengths), label="Length")
    plt.legend()
    plt.title("Infection over time")
    plt.xlabel("Time")
    plt.ylabel("[cm]")
    plt.show()