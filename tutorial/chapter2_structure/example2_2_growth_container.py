"""small example in a container"""

import numpy as np

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

# This example loops through a number of SDFs (cylinder etc.) and saves/plots them together with the root system
rs = pb.Plant()

# Open plant and root parameter from a file
path = path = "../../modelparameter/structural/rootsystem/"
name = "Zea_mays_4_Leitner_2014"
rs.readParameters(path + name + ".xml")

# 0. creates a cylindrical soil core with top radius 5 cm, bot radius 5 cm, height 50 cm, not square but circular
soilcore = pb.SDF_PlantContainer(5, 5, 40, False)  # |\label{l2_2_1:cylinder}|

# 1. Creates a square rhizotron r*r, with height h, rotated around the x-axis
r, h, alpha = 20, 4, 45
rhizotron = pb.SDF_PlantContainer(r, r, h, True)  # |\label{l2_2_1:rectangle}|
posA = pb.Vector3d(0, r, -h / 2)  # origin before rotation #|\label{l2_2_1:rotationbegin}|
A = pb.Matrix3d.rotX(alpha / 180.0 * np.pi)
posA = A.times(posA)  # origin after rotation
rotatedRhizotron = pb.SDF_RotateTranslate(rhizotron, alpha, 0, posA.times(-1))  # |\label{l2_2_1:rotationend}|

# 2. A split pot experiment
topBox = pb.SDF_PlantBox(22, 20, 5)  # |\label{l2_2_1:splitboxbegin}|
sideBox = pb.SDF_PlantBox(10, 20, 35)
left = pb.SDF_RotateTranslate(sideBox, pb.Vector3d(-6, 0, -5))
right = pb.SDF_RotateTranslate(sideBox, pb.Vector3d(6, 0, -5))
box_ = []
box_.append(topBox)
box_.append(left)
box_.append(right)
splitBox = pb.SDF_Union(box_)  # |\label{l2_2_1:splitboxend}|

# 3. Rhizotubes as obstacles
box = pb.SDF_PlantBox(96, 126, 130)  # box #|\label{l2_2_1:rhizobox}|
rhizotube = pb.SDF_PlantContainer(6.4, 6.4, 96, False)  # a single rhizotube #|\label{l2_2_1:rhizotube}|
rhizoX = pb.SDF_RotateTranslate(rhizotube, 90, pb.SDF_Axis.yaxis, pb.Vector3d(96 / 2, 0, 0))  # |\label{l2_2_1:rotatetube}|

rhizotubes_ = []  # |\label{l2_2_1:tubesmultbegin}|
y_ = (-30, -18, -6, 6, 18, 30)
z_ = (-10, -20, -40, -60, -80, -120)
tube = []
for i, y in enumerate(y_):
    v = pb.Vector3d(0, y, z_[i])
    tube.append(pb.SDF_RotateTranslate(rhizoX, v))
    rhizotubes_.append(tube[i])

rhizotubes = pb.SDF_Union(rhizotubes_)  # |\label{l2_2_1:tubesmultend}|
rhizoTube = pb.SDF_Difference(box, rhizotubes)  # |\label{l2_2_1:boxminustubes}|

#
number = 0
containers = [soilcore, rotatedRhizotron, splitBox, rhizoTube]
for container in containers:
    rs.setGeometry(container)
    rs.initialize()
    rs.simulate(45)  # days

    # Plot, using vtk
    vp.plot_roots_and_container(rs, container)

    # Export results (as vtp)
    rs.write("results/root_system" + str(number) + ".vtp")
    vp.write_container(container, "results/container_" + str(number) + ".vtp")  # you can pass , resolution=200 to increase resulting mesh quality, default is 100
    number += 1
