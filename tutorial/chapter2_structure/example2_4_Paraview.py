"""increase axial resolution (e.g. for animation)"""

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

plant = pb.Plant()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
plant.readParameters(path + name + ".xml")

# Modify axial resolution
for p in plant.getOrganRandomParameter(pb.root):
    p.dx = 0.1  # adjust resolution

# Simulate
plant.initialize()
plant.simulate(20)  # days

# Export results as segments
ana = pb.SegmentAnalyser(plant)
ana.write("results/animation.vtp")

# Export geometry as Paraview Python script
box = pb.SDF_PlantBox(15, 10, 40)
vp.write_container(box, "results/container.vtp")
