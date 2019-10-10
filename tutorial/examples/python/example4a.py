"""everything from scratch (without parameter files)"""
import sys
sys.path.append("../../..")
import plantbox as pb

import math

rs = pb.RootSystem()

# Root type parameter
p0 = pb.RootRandomParameter(rs)  # with default values,
p1 = pb.RootRandomParameter(rs)  # all standard deviations are 0

p0.name = "taproot"
p0.subType = 1
p0.lb = 1
p0.la = 10
p0.nob = 20
p0.ln = 89. / 19.
p0.theta = 30. / 180.*math.pi
p0.r = 1
p0.dx = 0.5
p0.successor = [2]  # add successors
p0.successorP = [1]
p0.tropismT = pb.TropismType.gravi
p0.tropismN = 1.
p0.tropismS = 0.2

p1.name = "lateral"
p1.subType = 2
p1.la = 25
p1.las = 10  # add standard deviation
p1.ln = 0
p1.r = 2
p1.dx = 0.1
p1.tropismS = 0.3

rs.setOrganRandomParameter(p0)
rs.setOrganRandomParameter(p1)

# Root system parameter (neglecting shoot borne)

rsp = pb.SeedRandomParameter(rs)
rsp.seedPos = pb.Vector3d(0., 0., -3.)
rsp.maxB = 100
rsp.firstB = 10.
rsp.delayB = 3.
rs.setRootSystemParameter(rsp)

rs.initialize()
rs.simulate(40, False)

rs.write("../results/example_3c.vtp")

print("done.")
