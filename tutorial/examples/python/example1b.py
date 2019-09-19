"""small example in a container"""
import py_rootbox as rb

rootsystem = rb.RootSystem()

# Open plant and root parameter from a file
name = "Anagallis_femina_Leitner_2010"
rootsystem.readParameters("modelparameter/" + name + ".xml")

# Create and set geometry

# 1.creates a cylindrical soil core with top radius 5 cm, bot radius 5 cm, height 50 cm, not square but circular
soilcore = rb.SDF_PlantContainer(5, 5, 40, False)

# 2. creates a square 27*27 cm containter with height 1.4 cm
rhizotron = rb.SDF_PlantBox(1.4, 27, 27)

# Pick 1, or 2
rootsystem.setGeometry(soilcore)  # soilcore, or rhizotron

# Initialize
rootsystem.initialize()

# Simulate
rootsystem.simulate(60)  # days

# Export final result (as vtp)
rootsystem.write("../results/example_1b.vtp")

# Export container geometry as Paraview Python script
rootsystem.write("../results/example_1b.py")

print("done.")
