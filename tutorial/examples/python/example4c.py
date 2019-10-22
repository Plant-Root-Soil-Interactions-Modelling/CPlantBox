"""increase axial resolution (e.g. for animation)"""
import sys
sys.path.append("../../..")
import plantbox as pb

rs = pb.RootSystem()

path = "../../../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rs.readParameters(path + name + ".xml")

# Set Geometry
soilcore = pb.SDF_PlantContainer(5, 5, 40, False)
rs.setGeometry(soilcore)

# Modify axial resolution
for p in rs.getRootRandomParameter():
    p.dx = 0.1  # adjust resolution

# Simulate
rs.initialize()
rs.simulate(60)  # days

# Export results as segments
ana = pb.SegmentAnalyser(rs)
ana.write("results/example_4c.vtp")

# Export container geometry as Paraview Python script
rs.write("results/example_4c.py")
