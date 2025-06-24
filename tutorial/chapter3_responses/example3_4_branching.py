"""scales branching probability"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np

plant = pb.Plant()
path = "../../modelparameter/structural/rootsystem/"
name = "Glycine_max_Moraes2020"
plant.readParameters(path + name + ".xml")

box = pb.SDF_PlantBox(10, 10, 30)  # nutrient rich patch  # |\label{l34:patch}|
patch = pb.SDF_RotateTranslate(box, pb.Vector3d(-5, 0., -10))

max_ = 1.  # maximal
min_ = 0.02  # minimal
slope = 1.  # [cm] linear gradient between min and max
soilprop = pb.SoilLookUpSDF(patch, max_, min_, slope)  # |\label{l34:rate}|

p = plant.getOrganRandomParameter(pb.root, 2)
p.ln = p.ln / 5  # increase overall branching density for subType 2 # |\label{l34:increase}|
p = plant.getOrganRandomParameter(pb.root, 3)
p.f_sbp = soilprop  # set branching probability for subType 3 # |\label{l34:prob}|

plant.initialize()  # |\label{l34:loop}|
simtime = 15.
dt = 1.
for i in range(0, round(simtime / dt)):
    plant.simulate(dt)  # |\label{l34:loop_end}|

l = plant.getSummed("length")  # |\label{l34:analysis}|
ana = pb.SegmentAnalyser(plant)
ana.crop(patch)
l_in = ana.getSummed("length")
ana = pb.SegmentAnalyser(plant)
ana.crop(pb.SDF_Complement(patch))
l_out = ana.getSummed("length")
print('\nRoot length within patch {:g} cm, {:g}%'.format(l_in, 100 * l_in / l))
print('Root length outside patch {:g} cm, {:g}% \n'.format(l_out, 100 * l_out / l))  # |\label{l34:analysis_end}|

plant.write("results/example_3_4a.vtp")
vp.plot_roots_and_container(plant, patch)
