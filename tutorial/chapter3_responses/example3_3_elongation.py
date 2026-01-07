"""scales the root elongation rate"""
<<<<<<< HEAD
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import numpy as np
=======

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
>>>>>>> origin/master

rs = pb.Plant()
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
<<<<<<< HEAD
maxS = 1.  # maximal
minS = 0.05  # minimal
slope = 1.  # [cm] linear gradient between min and max
=======
maxS = 1.0  # maximal
minS = 0.05  # minimal
slope = 1.0  # [cm] linear gradient between min and max
>>>>>>> origin/master
leftC = pb.SDF_Complement(left)
soilprop = pb.SoilLookUpSDF(leftC, maxS, minS, slope)

# Manually set scaling function and tropism parameters
for organ_type in [pb.root, pb.stem, pb.stem]:
<<<<<<< HEAD
    sigma = [0.4, 1., 1., 1., 1. ] * 2
=======
    sigma = [0.4, 1.0, 1.0, 1.0, 1.0] * 2
>>>>>>> origin/master
    for p in rs.getOrganRandomParameter(organ_type):
        p.dx = 0.25  # adjust resolutionx
        p.tropismS = sigma[p.subType - 1]
        p.f_se = soilprop  # 1. Scale elongation

# simulation
rs.initialize()
<<<<<<< HEAD
simtime = 60.
dt = 1.
=======
simtime = 60.0
dt = 1.0
>>>>>>> origin/master
for i in range(0, round(simtime / dt)):
    # in a dynamic setting change soilprop here
    rs.simulate(dt, False)

# analysis
l = rs.getSummed("length")
al = pb.SegmentAnalyser(rs)
al.crop(left)
ll = al.getSummed("length")
ar = pb.SegmentAnalyser(rs)
ar.crop(right)
lr = ar.getSummed("length")
<<<<<<< HEAD
print('\nLeft  compartment total root length {:g} cm, {:g}%'.format(ll, 100 * ll / l))
print('\nRight compartment total root length {:g} cm, {:g}% \n'.format(lr, 100 * lr / l))
=======
print(f"\nLeft  compartment total root length {ll:g} cm, {100 * ll / l:g}%")
print(f"\nRight compartment total root length {lr:g} cm, {100 * lr / l:g}% \n")
>>>>>>> origin/master

# write results
rs.write("results/example_5a.py")  # compartment geometry
rs.write("results/example_5a.vtp")  # root system

# plot, using vtk
vp.plot_roots(rs, "rootLength")  # press 'y'
