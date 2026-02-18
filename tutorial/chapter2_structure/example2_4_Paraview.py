"""increase axial resolution (e.g. for animation)"""
import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

plant = pb.Plant()
path = "../../modelparameter/structural/rootsystem/"
filename = "Anagallis_femina_Leitner_2010"
plant.readParameters(path + filename + ".xml")

# Modify axial resolution
for p in plant.getOrganRandomParameter(pb.root): # |\label{2_4:axres_start}|
    p.dx = 0.1  # adjust resolution              # |\label{2_4:axres_end}|

# Simulate
plant.initialize()                              
plant.simulate(20)  # days                      

# Export results as segments
ana = pb.SegmentAnalyser(plant)                 # |\label{2_4:seg_export_start}|
ana.write("results/animation.vtp")              # |\label{2_4:seg_export_end}|

# Export geometry as Paraview Python script
box = pb.SDF_PlantBox(15, 10, 40)                   # |\label{2_4:container_start}|
vp.write_container(box, "results/container.vtp")    # |\label{2_4:container_end}|
