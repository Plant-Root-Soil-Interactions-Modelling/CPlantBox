"""increase axial resolution (e.g. for animation)"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

plant = pb.Plant()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
plant.readParameters(path + name + ".xml")

# Modify axial resolution
for p in plant.getOrganRandomParameter(pb.root):
    p.dx = 0.1  # adjust resolution

# Simulate
plant.initialize()
plant.simulate(60)  # days

# Export results as segments
ana = pb.SegmentAnalyser(plant)
ana.write("results/animation.vtp")

ana.mapPeriodic(15, 10)
ana.write("results/animation_periodic.vtp")

# Export geometry as Paraview Python script
box = pb.SDF_PlantBox(15, 10, 50)
plant.setGeometry(box)
plant.write("results/animation_periodic.py")

# Plot final (periodic) image, using vtk
vp.plot_roots(ana, "creationTime")
