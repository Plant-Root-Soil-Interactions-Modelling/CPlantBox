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
p0.a = 0.2  # [cm] radius
p0.subType = 1 # [-] index starts at 1
p0.lb = 1  #[cm] basal zone
p0.la = 10 # [cm] apical zone 
p0.lmax = 100 # [cm] maximal root length, number of laterals = round((lmax-lb-la)/ln) + 1
p0.ln = 89. / 19. # [cm] inter-lateral distance 
p0.theta = 30. / 180.*math.pi # [rad]
p0.r = 1  # [cm/day] initial growth rate
p0.dx = 0.5 # [cm] axial resolution
p0.successor = [2]  # add successors
p0.successorP = [1] # probability that successor emerges
p0.tropismT = pb.TropismType.gravi #  
p0.tropismN = 1. # [-] strength of tropism
p0.tropismS = 0.2 # [rad/cm] maximal bending

p1.name = "lateral"
p1.a = 0.1 # [cm] radius
p1.subType = 2 # [1] index starts at 1
p1.la = 25 # # [cm] apical zone
p1.las = 10  # [cm] standard deviation of the apical zone
p1.ln = 0 # [cm] inter-lateral distance 
p1.r = 2  # initial growth rate
p1.dx = 0.1 # [cm] axial resolution
p1.tropismS = 0.3 # [rad/cm] maximal bending

rs.setOrganRandomParameter(p0)
rs.setOrganRandomParameter(p1)

# Seed random parameter (neglecting shoot borne)
srp = pb.SeedRandomParameter(rs) # with default values
srp.seedPos = pb.Vector3d(0., 0., -3.) # [cm] seed position
srp.maxB = 100 # [-] number of basal roots
srp.firstB = 10. # [day] first emergence of a basal root
srp.delayB = 3. # [day] delay between the emergence of basal roots
rs.setRootSystemParameter(srp)

rs.initialize(1, 1)  # basalType, shootborneType
rs.simulate(40) #  [day]
rs.write("../results/example_2a.vtp")

# Some output
print()
print("subType", rs.getParameter("subType"))
print("     la", rs.getParameter("la"))
print("la_mean", rs.getParameter("la_mean"))
print("radius", rs.getParameter("radius"))
print()
print("length", rs.getParameter("length"))
print("   age", rs.getParameter("age"))