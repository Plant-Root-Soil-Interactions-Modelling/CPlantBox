"""scales branching probability"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np

plant = pb.Plant()
path = "../../modelparameter/structural/rootsystem/"
name = "Glycine_max_Moraes2020"
plant.readParameters(path + name + ".xml")

box = pb.SDF_PlantBox(10, 10, 50)  # nutrient rich patch
patch = pb.SDF_RotateTranslate(box, pb.Vector3d(-4.99, 0, 0))

max_ = 1.  # maximal
min_ = 0.002  # minimal
slope = 1.  # [cm] linear gradient between min and max
soilprop = pb.SoilLookUpSDF(patch, max_, min_, slope)

# 3. Scale branching probability
p = plant.getOrganRandomParameter(pb.root, 2)
p.ln = p.ln / 5  # increase overall branching density
p = plant.getOrganRandomParameter(pb.root, 3)
p.f_sbp = soilprop

# Simulation
plant.initialize()
simtime = 15.
dt = 1.
for i in range(0, round(simtime / dt)):
    plant.simulate(dt)

# Analysis
l = plant.getSummed("length")
ana = pb.SegmentAnalyser(plant)
ana.crop(patch)
l_in = ana.getSummed("length")
ana = pb.SegmentAnalyser(plant)
ana.crop(pb.SDF_Complement(patch))
l_out = ana.getSummed("length")
print('\nLeft  compartment total root length {:g} cm, {:g}%'.format(l_in, 100 * l_in / l))
print('Right compartment total root length {:g} cm, {:g}% \n'.format(l_out, 100 * l_out / l))

plant.write("results/example_3_5.vtp")  # write results
vp.plot_roots_and_container(plant, patch)

