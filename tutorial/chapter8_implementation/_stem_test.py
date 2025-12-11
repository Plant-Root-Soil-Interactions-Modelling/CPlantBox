import numpy as np

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

plant = pb.Plant()

# test = pb.OrganRandomParameter(plant)
# print(test)

path = "../../modelparameter/structural/plant/"  # Open plant and root parameter from a file
name = "fspm2023"
plant.readParameters(path + name + ".xml")

srp = plant.getOrganRandomParameter(pb.seed)
rrp = plant.getOrganRandomParameter(pb.root)
strp = plant.getOrganRandomParameter(pb.stem)
lrp = plant.getOrganRandomParameter(pb.leaf)

# print(srp[0])
srp[0].firstSB = 1.e6
srp[0].firstB = 1.e6
srp[0].firstTil = 1.e6

# print(rrp[1])
# print(rrp[1].successor)
rrp[1].successor = [[]]
rrp[1].successorP = [[]]

# for i, s in enumerate(strp[1:]):
#     print("***", i)
#     print(s)
#     print()

# strp[0].successorOT = [[]]
# strp[0].successor = [[]]
# strp[0].successorP = [[]]
# strp[0].successorNo = []

print(strp[1])
strp[1].ln = 4  # BUG? 4
strp[1].lmax = 12  # BUG? lmax = 10
strp[1].RotBeta = 0
strp[1].BetaDev = 0
strp[1].InitBeta = 0

# strp[1].successorOT = [[]]
# strp[1].successor = [[]]
# strp[1].successorP = [[]]
# strp[1].successorNo = [1, 1]

strp[2].theta = np.pi / 2
# strp[2].RotBeta = 0
# strp[2].BetaDev = 0
# strp[2].InitBeta = 0
# strp[2].successorOT = [[]]
# strp[2].successor = [[]]
# strp[2].successorP = [[]]
# strp[2].successorP = []

lrp[1].dx = 0.05

plant.writeParameters(name + "_modified.xml")

plant.initialize()  # Initialize |\label{l13:initialize}|
simtime = 80  # days
plant.simulate(simtime)  # Simulate|\label{l13:simulate}|

ana = pb.SegmentAnalyser(plant)
vp.plot_plant(ana, "subType")  # e.g. organType, subType, age |\label{l13:plot_plant}|
