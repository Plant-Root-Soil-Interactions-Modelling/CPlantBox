"""increase axial resolution (e.g. for animation)"""
import sys
sys.path.append("../../..")
import plantbox as pb

rs = pb.RootSystem()
path = "../../../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rs.readParameters(path + name + ".xml")

# Modify axial resolution
for p in rs.getRootRandomParameter():
    p.dx = 0.1  # adjust resolution

# Simulate
rs.initialize()
rs.simulate(60)  # days

# Export results as segments
ana = pb.SegmentAnalyser(rs)
ana.write("results/example_4c.vtp")

ana.mapPeriodic(20, 15)
ana.write("results/example_4c_periodic.vtp")

# Export geometry as Paraview Python script
box = pb.SDF_PlantBox(20, 15, 35)
rs.setGeometry(box)
rs.write("results/example_4c_periodic.py")

