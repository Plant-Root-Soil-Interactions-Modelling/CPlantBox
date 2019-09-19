""""more complex geometries"""
import py_rootbox as rb
import math

rs = rb.RootSystem()

# Open plant and root parameter from a file
name = "Zea_mays_4_Leitner_2014"
rs.readParameters("modelparameter/" + name + ".xml")

# 1. Creates a square rhizotron r*r, with height h, rotated around the x-axis
r, h, alpha = 20, 4, 45
rhizotron2 = rb.SDF_PlantContainer(r, r, h, True)
posA = rb.Vector3d(0, r, -h / 2)  # origin before rotation
A = rb.Matrix3d.rotX(alpha / 180.*math.pi)
posA = A.times(posA)  # origin after rotation
rotatedRhizotron = rb.SDF_RotateTranslate(rhizotron2, alpha, 0, posA.times(-1))

# 2. A split pot experiment
topBox = rb.SDF_PlantBox(22, 20, 5)
sideBox = rb.SDF_PlantBox(10, 20, 35)
left = rb.SDF_RotateTranslate(sideBox, rb.Vector3d(-6, 0, -5))
right = rb.SDF_RotateTranslate(sideBox, rb.Vector3d(6, 0, -5))
box_ = rb.std_vector_SDF_()
box_.append(topBox)
box_.append(left)
box_.append(right)
splitBox = rb.SDF_Union(box_)

# 3. Rhizotubes as obstacles
box = rb.SDF_PlantBox(96, 126, 130)  # box
rhizotube = rb.SDF_PlantContainer(6.4, 6.4, 96, False)  # a single rhizotube
rhizoX = rb.SDF_RotateTranslate(rhizotube, 90, rb.SDF_Axis.yaxis, rb.Vector3d(96 / 2, 0, 0))

rhizotubes_ = rb.std_vector_SDF_()
y_ = (-30, -18, -6, 6, 18, 30)
z_ = (-10, -20, -40, -60, -80, -120)
tube = []
for i in range(0, len(y_)):
    v = rb.Vector3d(0, y_[i], z_[i])
    tube.append(rb.SDF_RotateTranslate(rhizoX, v))
    rhizotubes_.append(tube[i])

rhizotubes = rb.SDF_Union(rhizotubes_)
rhizoTube = rb.SDF_Difference(box, rhizotubes)

# Set geometry: rotatedRhizotron, splitBox, or rhizoTube
rs.setGeometry(rhizoTube)

# Simulate
rs.initialize()
rs.simulate(90)  # days

# Export results (as vtp)
rs.write("../results/example_2a.vtp")

# Export container geometry as Paraview Python script
rs.write("../results/example_2a.py")

print("done.")

