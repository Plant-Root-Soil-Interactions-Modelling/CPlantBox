"""everything from scratch (without parameter files)"""
import sys
sys.path.append("../../..")
import plantbox as pb

import math

rs = pb.RootSystem()

# Root random parameter
p0 = pb.RootRandomParameter(rs)  # with default values,
p1 = pb.RootRandomParameter(rs)  # all standard deviations are 0

p0.name = "taproot"
p0.a = 0.2  # cm radius
p0.subType = 1
p0.lb = 1
p0.la = 10
p0.nob = 20
p0.ln = 89. / 19.
p0.theta = 30. / 180.*math.pi
p0.r = 1  # initial growth rate
p0.dx = 0.5
p0.successor = [2]  # add successors
p0.successorP = [1]
p0.tropismT = pb.TropismType.gravi
p0.tropismN = 1.
p0.tropismS = 0.2

p1.name = "lateral"
p1.a = 0.1  # cm radius
p1.subType = 2
p1.la = 25
p1.las = 10  # add standard deviation
p1.ln = 0
p1.r = 2  # initial growth rate
p1.dx = 0.1
p1.tropismS = 0.3

rs.setOrganRandomParameter(p0)
rs.setOrganRandomParameter(p1)

# Seed random parameter (neglecting shoot borne)
srp = pb.SeedRandomParameter(rs)
srp.seedPos = pb.Vector3d(0., 0., -3.)
srp.maxB = 100
srp.firstB = 10.
srp.delayB = 3.
rs.setRootSystemParameter(srp)

rs.initialize(1, 1)  # basalType, shootborneType
rs.simulate(40)
rs.write("../results/example_4a.vtp")

# Some output
print()
print(rs.getParameter("subType"))
print(rs.getParameter("la"))
print(rs.getParameter("la_mean"))
print(rs.getParameter("radius"))
print()
print(rs.getParameter("length"))
print(rs.getParameter("age"))

