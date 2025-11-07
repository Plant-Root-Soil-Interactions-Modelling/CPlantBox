import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import math
import visualisation.vtk_plot as vp
import numpy as np

# make choice of treatment here

treatment = 1
animation = False

treatments = ["WS", "W", "RMCS", "RMC"]

names = ["TomatoJohanna_WildType_Silicone", "TomatoJohanna_WildType", "TomatoJohanna_RMC_Silicone", "TomatoJohanna_RMC"]

# total root lengths for each sample in mm after 108 days T5 to Q9 in Excel sheet under "TotalRootLength"
WSLengths = [76145.88048, 133110.4718, 96806.52112, 94722.62373, 60687.46504]
WLengths = [60645.95887, 70667.46043, 58991.74297, 88779.9162, 86080.82346]
RMCSLengths = [125165.1072, 134637.2097, 160151.4099, 143844.2341, 151240.8598]
RMCLengths = [97828.20896, 99849.61423, 108539.5156, 113054.6591, 101631.7984]

# average diameter per sample in mm after 108 days Q4 to T8 in Excel sheet under "Average diameter"

WSDia = [0.313235545, 0.294747188, 0.358293929, 0.351905364, 0.329237833]
WDia = [0.339339727, 0.342339, 0.342627182, 0.336897846, 0.324256143]
RMCSDia = [0.36291875, 0.295492267, 0.332729188, 0.313729263, 0.321435375]
RMCDia = [0.277527526, 0.352695111, 0.276951933, 0.295755846, 0.274213867]

expWLength = sum(WLengths)/len(WLengths)/10 # cm
expRMCLength = sum(RMCLengths)/len(RMCLengths)/10 # cm

expWDia = sum(WDia)/len(WDia)
expWRadius = expWDia / 20
expRMCDia = sum(RMCDia)/len(RMCDia)
expRMCRadius = expRMCDia / 20

expWLengthSD = math.sqrt(sum((x/10 - expWLength)**2 for x in WLengths)/len(WLengths))
expRMCLengthSD = math.sqrt(sum((x/10 - expRMCLength)**2 for x in RMCLengths)/len(RMCLengths))

expWDiaSD = math.sqrt(sum((x - expWDia)**2 for x in WDia)/len(WDia))
expWRadiusSD = math.sqrt(sum((x/20 - expWRadius)**2 for x in WDia)/len(WDia))
expRMCDiaSD = math.sqrt(sum((x - expRMCDia)**2 for x in RMCDia)/len(RMCDia))
expRMCRadiusSD = math.sqrt(sum((x/20 - expRMCRadius)**2 for x in RMCDia)/len(RMCDia))

# Expected values for comparison

WT = [expWRadius, expWRadiusSD, expWLength, expWLengthSD]
RMC = [expRMCRadius, expRMCRadiusSD, expRMCLength, expRMCLengthSD]

print("Expected values (mean radius cm, radius SD cm, mean length cm, length SD cm):")
if treatment == 1:
    print("Wild Type:", WT)
if treatment == 3:
    print("RMC:", RMC)
mycp = pb.MycorrhizalPlant()

path = ""
name = names[treatment]

print("Using parameter file:", path + name + ".xml")
mycp.readParameters(path +name + ".xml", fromFile = True, verbose = True)


# Pot dimensions 11.7 cm x 11.7 cm x 30 cm
pot = pb.SDF_PlantContainer(11.7, 11.7, 30, False)


mycp.setGeometry(pot)
mycp.initialize(True)

totalduration = 108
potduration = 73
simtime = totalduration
fpd = 1 # one step per day
N = simtime * fpd
dt = simtime / N

for i in range(1,N):
        if i % 18 == 0:
            print('step',i, '/',N)
        mycp.simulate(dt, False)
        if animation:
            ana = pb.SegmentAnalyser(mycp)
            ana.addData("infection", mycp.getNodeInfections(2))
            ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
            ana.write(names[treatment] + "_animation_" + "{:04d}".format(i) + ".vtp", ["radius", "subType", "creationTime", "length", "infection", "infectionTime"])
    
print('done')

ana = pb.SegmentAnalyser(mycp)
rad = ana.getParameter("radius")
rad = sum(rad)/len(rad)
rl = ana.getSummed("length")
print("\n")
if treatment == 1:
    print("Simulated Wild Type:")
    print("Root length total (cm):", rl)
    print("Difference to expected lengths (cm):", rl - WT[2])
    print("Standard deviation of expected lengths (cm):", WT[3])
    print("Average root radius (cm):", rad)
    print("Difference to expected radius (cm):", rad - WT[0])
    print("Standard deviation of expected radius (cm):", WT[1])

if treatment == 3:
    print("Simulated AMF Resistant Type:")
    print("Root length total (cm):", rl)
    print("Difference to expected lengths (cm):", rl - RMC[2])
    print("Standard deviation of expected lengths (cm):", RMC[3])
    print("Average root radius (cm):", rad)
    print("Difference to expected radius (cm):", rad - RMC[0])
    print("Standard deviation of expected radius (cm):", RMC[1])


vp.plot_roots_and_container(mycp, pot)
ana.write(names[treatment] + ".vtp", ["radius", "subType", "creationTime","organType"])
vp.write_container(pot, "pot.vtp")


