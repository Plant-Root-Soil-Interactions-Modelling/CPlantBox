"""scales insertion angle"""
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
minS = 0.1  # minimal
slope = 1.  # [cm] linear gradient between min and max
leftC = pb.SDF_Complement(left)
soilprop = pb.SoilLookUpSDF(leftC, maxS, minS, slope)

# Manually set scaling function and tropism parameters
sigma = [0.4, 1., 1., 1., 1. ] * 2
for p in rs.getRootRandomParameter():
    if p.subType > 2:
        p.dx = 0.25  # adjust resolution
        p.f_sa = soilprop  # Scale insertion angle
        p.lmax = 2 * p.lmax  # make second order laterals longer

# simulation
rs.initialize()
simtime = 60.
dt = 1.
for i in range(0, round(simtime / dt)):
    rs.simulate(dt, False)

# analysis
al = pb.SegmentAnalyser(rs)
al.crop(left)
lm_theta = np.mean(al.getParameter("theta"))
ar = pb.SegmentAnalyser(rs)
ar.crop(right)
rm_theta = np.mean(ar.getParameter("theta"))
print('\nLeft  compartment mean insertion angle is {:g} degrees'.format(lm_theta))
print('\nRight compartment mean insertion angle is {:g} degrees\n'.format(rm_theta))

# write results
rs.write("results/example_5c.py")  # compartment geometry
rs.write("results/example_5c.vtp")  # root system

# plot, using vtk
vp.plot_roots(rs, "theta")  # press 'y'
