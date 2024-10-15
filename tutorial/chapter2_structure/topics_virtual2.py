""""more complex plant containers"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np

plant = pb.Plant()
path = path = "../../modelparameter/structural/rootsystem/"
name = "Zea_mays_4_Leitner_2014"
plant.readParameters(path + name + ".xml")

# 1. Creates a square rhizotron r*r, with height h, rotated around the x-axis
r, h, alpha = 20, 4, 45
rhizotron2 = pb.SDF_PlantContainer(r, r, h, True)
posA = pb.Vector3d(0, r, -h / 2)  # seed location bevore rotation
A = pb.Matrix3d.rotX(alpha / 180.*np.pi)
posA = A.times(posA)  # seed location after rotation
rotatedRhizotron = pb.SDF_RotateTranslate(rhizotron2, alpha, 0, posA.times(-1))

# 2. A split pot experiment
topBox = pb.SDF_PlantBox(22, 20, 5)
sideBox = pb.SDF_PlantBox(10, 20, 35)
left = pb.SDF_RotateTranslate(sideBox, pb.Vector3d(-6, 0, -5))
right = pb.SDF_RotateTranslate(sideBox, pb.Vector3d(6, 0, -5))
box_ = []
box_.append(topBox)
box_.append(left)
box_.append(right)
splitBox = pb.SDF_Union(box_)

# Simulate
plant.setGeometry(rotatedRhizotron)  # rotatedRhizotron, splitBox
plant.initialize()
plant.simulate(40)  # days

# Export and plot
plant.write("results/topics_virtual2.vtp")
plant.write("results/topics_virtual2.py")
vp.plot_roots(plant, "subType")
