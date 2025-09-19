"""small example in a cylindrical container or rhizotron"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

plant = pb.Plant()

path = path = "../../modelparameter/structural/rootsystem/"
name = "Zea_mays_4_Leitner_2014"
plant.readParameters(path + name + ".xml")

# 1.Creates a cylindrical container with top radius 5 cm, bot radius 5 cm, height 50 cm, not square but circular
soilcore = pb.SDF_PlantContainer(5, 5, 40, False)

# 2. Creates a square 27*27 cm containter with height 1.4 cm
rhizotron = pb.SDF_PlantBox(1.4, 27, 27)

# Pick 1, or 2
plant.setGeometry(rhizotron)  # soilcore, or rhizotron

# Initialize & Simulate
plant.initialize()
plant.simulate(40)  # days

# Export final result (as vtp)
plant.write("results/topics_virtual.vtp")
pb.SegmentAnalyser(plant).write("results/topics_virtual_seg.vtp")

# Export containr geometry as Paraview Python script
plant.write("results/topics_virtual.py")

# Plot, using vtk
vp.plot_plant(plant, "subType")
