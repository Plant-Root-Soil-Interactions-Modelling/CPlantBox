""" something basic"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np

plant = pb.Plant()
root = pb.RootRandomParameter(plant)
stem = pb.StemRandomParameter(plant)
leaf = pb.LeafRandomParameter(plant)
seed = pb.SeedRandomParameter(plant)

print(stem)
stem.subType = 1
stem.lmax = 10
stem.theta = 0.
leaf.subType = 1
root.subType = 1
root.lmax = 100

plant.setOrganRandomParameter(stem)
plant.setOrganRandomParameter(leaf)
plant.setOrganRandomParameter(seed)
plant.setOrganRandomParameter(root)

plant.initialize()
plant.simulate(20)
vp.plot_plant(plant, "organType")
