"""scales branching probability"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import numpy as np

rs = pb.RootSystem()
path = "../../modelparameter/structural/rootsystem/"
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
minS = 0.002  # minimal
slope = 1.  # [cm] linear gradient between min and max
leftC = pb.SDF_Complement(left)
soilprop = pb.SoilLookUpSDF(left, maxS, minS, slope)

# Manually set scaling function and tropism parameters
sigma = [0.4, 1., 1., 1., 1. ] * 2
for p in rs.getRootRandomParameter():
    p.dx = 0.25  # adjust resolution
    p.tropismS = sigma[p.subType - 1]

# 3. Scale branching probability
p = rs.getRootRandomParameter(2)
p.ln = p.ln / 5
p = rs.getRootRandomParameter(3)
p.f_sbp = soilprop

# simulation
rs.initialize()
simtime = 60.
dt = 1.
for i in range(0, round(simtime / dt)):
    rs.simulate(dt, True)

# analysis
l = rs.getSummed("length")
al = pb.SegmentAnalyser(rs)
al.crop(left)
ll = al.getSummed("length")
ar = pb.SegmentAnalyser(rs)
ar.crop(right)
lr = ar.getSummed("length")
print('\nLeft  compartment total root length {:g} cm, {:g}%'.format(ll, 100 * ll / l))
print('\nRight compartment total root length {:g} cm, {:g}% \n'.format(lr, 100 * lr / l))

# write results
rs.write("results/example_5d.py")  # compartment geometry
rs.write("results/example_5d.vtp")  # root system

# Plot, using vtk
vp.plot_roots(rs, "type")
