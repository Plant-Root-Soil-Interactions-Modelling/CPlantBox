"""root system length over time"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import rsml.rsml_reader as rsml

import numpy as np
import matplotlib.pyplot as plt

path = ""
name = "WurzelAtlasParasTomate"
name_new = "WurzelAtlasTomate_parameters_out"
name_Johanna_WildType = "TomatoJohanna_WildType"
name_Johanna_RMC = "TomatoJohanna_RMC"

animation = False
pot = False;



tomato_WAT = pb.MycorrhizalPlant()
tomato_J_WT = pb.MycorrhizalPlant()
tomato_J_RMC = pb.MycorrhizalPlant()

tomato_WAT.readParameters(path +name_new + ".xml", fromFile = True, verbose = True)
tomato_J_WT.readParameters(path +name_Johanna_WildType + ".xml", fromFile = True, verbose = True)
tomato_J_RMC.readParameters(path +name_Johanna_RMC + ".xml", fromFile = True, verbose = True)

if pot:
    # Pot dimensions 11.7 cm x 11.7 cm x 30 cm
    pot = pb.SDF_PlantContainer(11.7, 11.7, 30, False)
    tomato_WAT.setGeometry(pot)
    tomato_J_WT.setGeometry(pot)
    tomato_J_RMC.setGeometry(pot) 

# root_paras = pb.HyphaeRandomParameter(tomato_WAT)
# root_paras.hyphalEmergenceDensity = 0
# root_paras.highresolution = 0
# tomato_WAT.setOrganRandomParameter(root_paras)

hyphae_parameter = pb.HyphaeRandomParameter(tomato_WAT)
hyphae_parameter.subType = 1
hyphae_parameter.a = 0.01
hyphae_parameter.dx = 0.01
tomato_WAT.setOrganRandomParameter(hyphae_parameter)

tomato_WAT.initialize(True)
tomato_J_WT.initialize(True)
tomato_J_RMC.initialize(True)

# rootWATParas = pb.MycorrhizalRootRandomParameter(tomato_WAT)
# rootJWTParas = pb.MycorrhizalRootRandomParameter(tomato_J_WT)
# rootJRMCParas = pb.MycorrhizalRootRandomParameter(tomato_J_RMC)
# print(rootWATParas)
# print(rootJWTParas)
# print(rootJRMCParas)


# mycp = pb.MycorrhizalPlant()
# mycp.readParameters(path +name + ".xml", fromFile = True, verbose = True)

# # mycp.setGeometry(pot)
# mycp.initialize(True)
# mycp.writeParameters("WurzelAtlasTomate_parameters_out.xml")

totalduration = 108
potduration = 73
simtime = totalduration
fpd = 1 # one step per day
N = simtime * fpd
dt = simtime / N

for i in range(1,N):
    if i % 18 == 0:
        print('step',i, '/',N)
    tomato_WAT.simulate(dt, True)
    tomato_J_WT.simulate(dt, True)
    tomato_J_RMC.simulate(dt, True)
        # mycp.simulate(dt, False)
print('done')

ana_WAT = pb.SegmentAnalyser(tomato_WAT)
ana_J_WT = pb.SegmentAnalyser(tomato_J_WT)
ana_J_RMC = pb.SegmentAnalyser(tomato_J_RMC)

rad_WAT = ana_WAT.getParameter("radius")
rad_WAT = sum(rad_WAT)/len(rad_WAT)
rl_WAT = ana_WAT.getSummed("length")
rad_J_WT = ana_J_WT.getParameter("radius")
rad_J_WT = sum(rad_J_WT)/len(rad_J_WT)
rl_J_WT = ana_J_WT.getSummed("length")
rad_J_RMC = ana_J_RMC.getParameter("radius")
rad_J_RMC = sum(rad_J_RMC)/len(rad_J_RMC)
rl_J_RMC = ana_J_RMC.getSummed("length")    

print("\n")

print("Simlated average root radius and total root length after", simtime, "days:", "\n")
print("WurzelAtlas Paras Tomato:", "\n","  Average root radius (cm):", rad_WAT, "\n","  Total root length (cm):", rl_WAT, "\n")
print("Johanna Wild Type Tomato:", "\n","  Average root radius (cm):", rad_J_WT, "\n","  Total root length (cm):", rl_J_WT, "\n")
print("Johanna RMC Tomato:", "\n","  Average root radius (cm):", rad_J_RMC, "\n","  Total root length (cm):", rl_J_RMC, "\n")  


ana_WAT.write(name_new + ".vtp", ["radius", "subType", "creationTime","organType"])
ana_J_WT.write(name_Johanna_WildType + ".vtp", ["radius", "subType", "creationTime","organType"])
ana_J_RMC.write(name_Johanna_RMC + ".vtp", ["radius", "subType", "creationTime","organType"])
# ana = pb.SegmentAnalyser(mycp)
# rad = ana.getParameter("radius")
# rad = sum(rad)/len(rad)
# rl = ana.getSummed("length")
# print("\n")
# print("Simulated Root Atlas Type:")
# print("Root length total (cm):", rl)
# # print("Difference to expected lengths (cm):", rl - WT[2])
# # print("Standard deviation of expected lengths (cm):", WT[3])
# print("Average root radius (cm):", rad)
# # print("Difference to expected radius (cm):", rad - WT[0])
# # print("Standard deviation of expected radius (cm):", WT[1])

# vp.plot_roots(mycp)
# ana.write(name + ".vtp", ["radius", "subType", "creationTime","organType"])
# # vp.write_container(pot, "pot.vtp")

