""" something basic"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np

plant = pb.Plant()
seed = pb.SeedRandomParameter(plant)
root = pb.RootRandomParameter(plant)
stem = pb.StemRandomParameter(plant)
leaf = pb.LeafRandomParameter(plant)

stem.subType = 1
stem.lmax = 10
stem.la, stem.lb, stem.ln = 1, 1, 1
stem.theta = 0.

stem.successorOT = [[pb.leaf]]
stem.successorST = [[1]]
stem.successorP = [[1]]

leaf.subType = 1
leaf.lmax = 10

leaf.shapeType = 2
leaf.areaMax = 10.
leaf.leafGeometryPhi = [0., np.pi / 2, np.pi, 3 * n.pi / 2, 2.*np.pi]
leaf.leafGeometryX = [1., 1., 1., 1., 1.]
leaf.parametrisationType = 1
print(leaf)

root.subType = 1
root.lmax = 100

plant.setOrganRandomParameter(seed)
plant.setOrganRandomParameter(stem)
plant.setOrganRandomParameter(leaf)
plant.setOrganRandomParameter(root)

plant.initialize()
plant.simulate(20)
vp.plot_plant(plant, "organType")
