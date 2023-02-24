""""more complex geometries"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import numpy as np

rs = pb.RootSystem()

# Open plant and root parameter from a file
path = path = "../../modelparameter/structural/rootsystem/"
name = "Zea_mays_4_Leitner_2014"
rs.readParameters(path + name + ".xml")

# 1. Creates a square rhizotron r*r, with height h, rotated around the x-axis
r, h, alpha = 20, 4, 45
rhizotron2 = pb.SDF_PlantContainer(r, r, h, True)
posA = pb.Vector3d(0, r, -h / 2)  # origin before rotation
A = pb.Matrix3d.rotX(alpha / 180.*np.pi)
posA = A.times(posA)  # origin after rotation
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

# 3. Rhizotubes as obstacles
box = pb.SDF_PlantBox(96, 126, 130)  # box
rhizotube = pb.SDF_PlantContainer(6.4, 6.4, 96, False)  # a single rhizotube
rhizoX = pb.SDF_RotateTranslate(rhizotube, 90, pb.SDF_Axis.yaxis, pb.Vector3d(96 / 2, 0, 0))

rhizotubes_ = []
y_ = (-30, -18, -6, 6, 18, 30)
z_ = (-10, -20, -40, -60, -80, -120)
tube = []
for i in range(0, len(y_)):
    v = pb.Vector3d(0, y_[i], z_[i])
    tube.append(pb.SDF_RotateTranslate(rhizoX, v))
    rhizotubes_.append(tube[i])

rhizotubes = pb.SDF_Union(rhizotubes_)
rhizoTube = pb.SDF_Difference(box, rhizotubes)

# Set geometry: rotatedRhizotron, splitBox, or rhizoTube
rs.setGeometry(splitBox)

# Simulate
rs.initialize()
rs.simulate(90)  # days

# Export results (as vtp)
rs.write("results/example_1c.vtp")

# Export container geometry as Paraview Python script
rs.write("results/example_1c.py")

# Plot, using vtk
vp.plot_roots(rs, "type")
