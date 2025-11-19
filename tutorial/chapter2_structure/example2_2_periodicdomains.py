"""increase axial resolution (e.g. for animation)"""

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

rs = pb.Plant()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rs.readParameters(path + name + ".xml")

# Modify axial resolution
for p in rs.getOrganRandomParameter(pb.root):
    p.dx = 0.1  # adjust resolution

# Simulate
rs.initialize()
rs.simulate(60)  # days

# Export results as segments
ana = pb.SegmentAnalyser(rs)
ana.write("results/periodic.vtp")

ana.mapPeriodic(15, 10)  # |\label{l2_2_3:periodicity}|
ana.write("results/periodic.vtp")

# Export geometry as Paraview Python script
box = pb.SDF_PlantBox(15, 10, 35)
rs.setGeometry(box)
rs.write("results/periodic.py")

# Plot final (periodic) image, using vtk
vp.plot_roots(ana, "creationTime")
