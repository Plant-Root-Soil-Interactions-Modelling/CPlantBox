"""Three types of interaction, setting f_se, f_sa, f_sbp"""
import sys
sys.path.append("../../..")
import plantbox as pb
import math

rs = pb.RootSystem()
path = "../../../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rs.readParameters(path + name + ".xml")

# box with a left and a right compartment for analysis
sideBox = pb.SDF_PlantBox(10, 20, 50)
left = pb.SDF_RotateTranslate(sideBox, pb.Vector3d(-4.99, 0, 0))
right = pb.SDF_RotateTranslate(sideBox, pb.Vector3d(4.99, 0, 0))
leftright = pb.SDF_Union(left, right)
rs.setGeometry(leftright)

# left compartment has a minimum of 0.01, 1 elsewhere
maxS = 1.  # maximal
minS = 0.01  # minimal
slope = 1.  # [cm] linear gradient between min and max
leftC = pb.SDF_Complement(left)
soilprop = pb.SoilLookUpSDF(leftC, maxS, minS, slope)  # for root elongation
soilprop2 = pb.SoilLookUpSDF(left, 1., 0.002, slope)  # for branching

# Manually set scaling function and tropism parameters
sigma = [0.4, 1., 1., 1., 1. ] * 2
for p in rs.getRootRandomParameter():
    p.dx = 0.25  # adjust resolution
    p.tropismS = sigma[p.subType - 1]

    # 1. Scale elongation
#    p.f_se = soilprop

    # 2. Scale insertion angle
#    p.f_sa = soilprop

# 3. Scale branching probability
p = rs.getRootRandomParameter(2)
p.ln = p.ln / 5
# p.nob = p.nob * 5
p = rs.getRootRandomParameter(3)
p.f_sbp = soilprop2

# simulation
rs.initialize()
simtime = 120.
dt = 1.
N = 120 / dt
for i in range(0, round(N)):
    rs.simulate(dt, True)

# analyse
print()
print("Left compartment: ")
al = pb.SegmentAnalyser(rs)
al.crop(left)
ll = al.getSummed("length")
print('Total root length', ll, 'cm')
lmct = al.getSummed("creationTime") / al.getSummed("one")
print('Mean age', simtime - lmct, 'days')
lroots = al.getOrgans()
lm_theta = 0
for r in lroots:
    lm_theta += r.param().theta  # downcast
lm_theta /= len(lroots)
print('Mean insertion angle is ', lm_theta / math.pi * 180, 'degrees')
print()

print("Right compartment: ")
ar = pb.SegmentAnalyser(rs)
ar.crop(right)
lr = ar.getSummed("length")
print('Total root length', lr, 'cm')
rmct = ar.getSummed("creationTime") / ar.getSummed("one")
print('Mean age', simtime - rmct, 'days')
rroots = ar.getOrgans()
rm_theta = 0
for r in lroots:
    rm_theta += r.param().theta
rm_theta /= len(rroots)
print('Mean insertion angle is ', rm_theta / math.pi * 180, 'degrees')
print()

# write results
rs.write("results/example_6a.py")  # compartment geometry
rs.write("results/example_6a.vtp")  # root system

