"""small example in a container"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

rs = pb.Plant()

# Open plant and root parameter from a file
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rs.readParameters(path + name + ".xml")

# Create and set geometry

# 1.creates a cylindrical soil core with top radius 5 cm, bot radius 5 cm, height 50 cm, not square but circular
soilcore = pb.SDF_PlantContainer(5, 5, 40, False)

# 2. creates a square 27*27 cm containter with height 1.4 cm
rhizotron = pb.SDF_PlantBox(1.4, 27, 27)

# Pick 1, or 2
rs.setGeometry(rhizotron)  # soilcore, or rhizotron

# Initialize
rs.initialize()

# Simulate
rs.simulate(60)  # days

# Export final result (as vtp)
rs.write("results/example_1b.vtp")

# Export container geometry as Paraview Python script
# rs.write("results/example_1b.py")

# Plot, using vtk
vp.plot_roots(rs, "type")
