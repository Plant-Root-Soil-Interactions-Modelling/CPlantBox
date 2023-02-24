"""increase axial resolution (e.g. for animation)"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

rs = pb.RootSystem()
path = "../../modelparameter/structural/rootsystem/"
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
ana.write("results/example_3e.vtp")

ana.mapPeriodic(15, 10)
ana.write("results/example_3e_periodic.vtp")

# Export geometry as Paraview Python script
box = pb.SDF_PlantBox(15, 10, 35)
rs.setGeometry(box)
rs.write("results/example_3e_periodic.py")

# Plot final (periodic) image, using vtk
vp.plot_roots(ana, "creationTime")
