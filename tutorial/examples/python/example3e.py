"""increase axial resolution (e.g. for animation)"""
import py_rootbox as rb

rs = rb.RootSystem()

# Open plant and root parameter from a file
name = "Anagallis_femina_Leitner_2010"
rs.readParameters("modelparameter/" + name + ".xml")

# Set Geometry
soilcore = rb.SDF_PlantContainer(5, 5, 40, False)
rs.setGeometry(soilcore)

# Modify axial resolution
for p in rs.getRootTypeParameter():
    p.dx = 0.1  # adjust resolution

# Simulate
rs.initialize()
rs.simulate(60, True)  # days

# Export results as segments
rb.SegmentAnalyser(rs).write("../results/example_3e.vtp")

# Export container geometry as Paraview Python script
rs.write("../results/example_3e.py")

print("done.")
