""""more complex plant containers with obstacles"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np

plant = pb.Plant()
path = path = "../../modelparameter/structural/rootsystem/"
name = "Moraesetal_2020"
plant.readParameters(path + name + ".xml")

# Rhizotubes as obstacles
box = pb.SDF_PlantBox(96, 126, 130)  # box
rhizotube = pb.SDF_PlantContainer(3., 3., 96, False)  # a single rhizotube
rhizoX = pb.SDF_RotateTranslate(rhizotube, 90, pb.SDF_Axis.yaxis, pb.Vector3d(96 / 2, 0, 0))

rhizotubes_ = []
y_ = (-30, -18, -6, 6, 18, 30)  # cm
z_ = (-10, -20, -40, -60, -80, -120)  # cm
tube = []
for i in range(0, len(y_)):
    v = pb.Vector3d(0, y_[i], z_[i])
    tube.append(pb.SDF_RotateTranslate(rhizoX, v))
    rhizotubes_.append(tube[i])

rhizotubes = pb.SDF_Union(rhizotubes_)
rhizoTube = pb.SDF_Difference(box, rhizotubes)

# Simulate
plant.setGeometry(rhizoTube)
plant.initialize()
plant.simulate(90)  # days

# Export results
plant.write("results/topics_virtual3.vtp")
plant.write("results/topics_virtual3.py")
vp.plot_roots(plant, "age")

